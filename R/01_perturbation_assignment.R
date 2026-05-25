#!/usr/bin/env Rscript
# =============================================================================
# CIRI Pipeline - Step 01: Perturbation Assignment
# =============================================================================
# Assigns CRISPR guide identities to each cell from a 10x Genomics H5 matrix.
# This is a Docker-free refactor of perturbation_assignment/CIRI_unified.R.
#
# Usage:
#   Rscript 01_perturbation_assignment.R \
#     --dir          /path/to/data/ \
#     --matrix       filtered_feature_bc_matrix.h5 \
#     --guides       guides.csv \
#     --strategy     1 \
#     --threshold_a  -1 \
#     --threshold_i  -1
#
# Arguments:
#   --dir          Directory containing the H5 matrix and guides.csv
#   --matrix       H5 matrix filename (inside --dir)
#   --guides       Guide annotation CSV (inside --dir). Default: guides.csv
#   --strategy     1 = single variable guide, 2 = dual variable guide
#   --threshold_a  UMI threshold for CRISPRa. -1 = auto-detect via KDE.
#   --threshold_i  UMI threshold for CRISPRi. -1 = auto-detect via KDE.
#
# Output (written to --dir):
#   annotation_data.csv        – per-cell guide assignment
#   CIRI_long.csv              – long-format UMI table
#   *.pdf                      – QC plots
# =============================================================================

suppressMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(quantmod)
  library(pracma)
  library(hdf5r)
  library(Matrix)
  library(zoo)
  library(scales)
  library(Seurat)
})

# ---------------------------------------------------------------------------
# Argument Parsing
# ---------------------------------------------------------------------------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  # Defaults
  params <- list(
    dir         = NULL,
    matrix      = NULL,
    guides      = "guides.csv",
    strategy    = 1L,
    threshold_a = -1,
    threshold_i = -1
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    val <- if (i + 1 <= length(args)) args[i + 1] else NA
    switch(key,
      "--dir"         = { params$dir         <- val; i <- i + 2 },
      "--matrix"      = { params$matrix       <- val; i <- i + 2 },
      "--guides"      = { params$guides       <- val; i <- i + 2 },
      "--strategy"    = { params$strategy     <- as.integer(val); i <- i + 2 },
      "--threshold_a" = { params$threshold_a  <- as.numeric(val); i <- i + 2 },
      "--threshold_i" = { params$threshold_i  <- as.numeric(val); i <- i + 2 },
      { message(paste("Unknown argument:", key)); i <- i + 1 }
    )
  }

  if (is.null(params$dir))    stop("--dir is required")
  if (is.null(params$matrix)) stop("--matrix is required")
  if (!params$strategy %in% c(1L, 2L)) stop("--strategy must be 1 or 2")

  params$dir <- normalizePath(params$dir, mustWork = TRUE)
  params
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

#' Auto-detect UMI threshold using the first derivative of KDE.
#' Finds the valley between the noise peak and the signal peak.
auto_threshold <- function(umi_vec, label = "") {
  umi_vec <- umi_vec[umi_vec > 0]
  if (length(umi_vec) < 20) {
    warning(paste("Too few positive UMI values to auto-threshold for", label, ". Using 5."))
    return(5)
  }

  dens   <- density(umi_vec, bw = "SJ", n = 1024)
  dy     <- diff(dens$y)
  # Find where derivative crosses zero from negative to positive (local minimum)
  valley_idx <- which(dy[-length(dy)] < 0 & dy[-1] >= 0)

  if (length(valley_idx) == 0) {
    warning(paste("Could not detect valley for", label, ". Using median threshold."))
    return(median(umi_vec) / 2)
  }

  thresh <- dens$x[valley_idx[1] + 1]
  message(paste("Auto-threshold for", label, ":", round(thresh, 2)))
  thresh
}

#' Load the CRISPR capture matrix from an HDF5 file.
load_crispr_h5 <- function(h5_path) {
  f        <- H5File$new(h5_path, mode = "r")
  on.exit(f$close_all(), add = TRUE)

  feat     <- f[["matrix/features"]]
  all_names  <- feat[["name"]][]
  types      <- feat[["feature_type"]][]
  crispr_idx <- which(types == "CRISPR Guide Capture")

  if (length(crispr_idx) == 0) stop("No 'CRISPR Guide Capture' features found in H5 file.")

  data     <- f[["matrix/data"]][]
  indices  <- f[["matrix/indices"]][]
  indptr   <- f[["matrix/indptr"]][]
  shape    <- f[["matrix/shape"]][]
  barcodes <- f[["matrix/barcodes"]][]

  mat          <- new("dgCMatrix", x = as.numeric(data), i = indices, p = indptr, Dim = shape)
  rownames(mat) <- all_names
  colnames(mat) <- barcodes

  mat[crispr_idx, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# QC Plots
# ---------------------------------------------------------------------------

save_qc_plots <- function(CIRI_long, crispra, crispri, dir) {
  # Total UMIs per guide
  plot_df <- CIRI_long %>%
    group_by(feature, type) %>%
    summarise(total_umi = sum(umi), .groups = "drop") %>%
    arrange(type, feature) %>%
    mutate(feature = factor(feature, levels = unique(feature)))

  p <- ggplot(plot_df, aes(x = feature, y = total_umi, fill = type)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Total UMIs per guide", x = "Guide", y = "Total UMIs") +
    scale_y_continuous(labels = label_number(accuracy = 1)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(file.path(dir, "total_umixguide.pdf"), p)

  # UMI histogram (all guides)
  p <- ggplot(CIRI_long, aes(x = umi)) + geom_histogram(binwidth = 1) +
    labs(title = "UMI distribution per guide")
  ggsave(file.path(dir, "cellranger_UMIxguide_hist.pdf"), p)

  # UMI per cell
  umis_per_cell <- CIRI_long %>%
    group_by(cell_barcode) %>%
    summarise(total_umis = sum(umi), .groups = "drop")
  p <- ggplot(umis_per_cell, aes(x = total_umis)) + geom_histogram(binwidth = 1) +
    labs(title = "Total UMIs per cell")
  ggsave(file.path(dir, "cellranger_UMIxcell_hist.pdf"), p)

  # Per-guide histograms
  dir.create(file.path(dir, "single_guide_plots"), showWarnings = FALSE)
  for (guide in unique(CIRI_long$feature)) {
    df <- filter(CIRI_long, feature == guide)
    ggsave(file.path(dir, "single_guide_plots", paste0("hist_umi_", guide, ".pdf")),
           ggplot(df, aes(x = umi)) + geom_histogram(binwidth = 1) + labs(title = guide))
    ggsave(file.path(dir, "single_guide_plots", paste0("hist_pct_", guide, ".pdf")),
           ggplot(df, aes(x = percentage)) + geom_histogram(binwidth = 1) + labs(title = guide))
  }

  # Fixed vs variable scatter
  CRISPRai_plot <- bind_rows(crispra, crispri) %>%
    mutate(variable_umi = total_umi - fixed_umi)
  p <- ggplot(CRISPRai_plot, aes(x = variable_umi, y = fixed_umi, color = feature)) +
    geom_point(alpha = 0.1, size = 0.5) + theme_minimal() +
    labs(title = "Fixed vs Variable UMI per Cell",
         x = "Variable UMI", y = "Fixed UMI", color = "Feature")
  ggsave(file.path(dir, "CRISPRai_fixed_vs_variable_scatter.pdf"), p, width = 12, height = 10)
}

# ---------------------------------------------------------------------------
# Assignment Logic - Single Variable Guide (strategy = 1)
# ---------------------------------------------------------------------------

assign_single <- function(CIRI_long, a_genes, i_genes, thresh_a, thresh_i) {
  # CRISPRa assignment
  crispra_df <- CIRI_long %>%
    filter(type == "a") %>%
    mutate(is_fixed = feature %in% a_genes) %>%
    group_by(cell_barcode) %>%
    summarise(
      total_umi  = sum(umi),
      fixed_umi  = sum(umi[is_fixed]),
      percentage = (fixed_umi / total_umi) * 100,
      feature    = "CRISPRa",
      .groups    = "drop"
    ) %>%
    filter(total_umi > 0)

  crispri_df <- CIRI_long %>%
    filter(type == "i") %>%
    mutate(is_fixed = feature %in% i_genes) %>%
    group_by(cell_barcode) %>%
    summarise(
      total_umi  = sum(umi),
      fixed_umi  = sum(umi[is_fixed]),
      percentage = (fixed_umi / total_umi) * 100,
      feature    = "CRISPRi",
      .groups    = "drop"
    ) %>%
    filter(total_umi > 0)

  # Identify positive cells (fixed guide present above threshold)
  pos_a_cells <- crispra_df %>% filter(fixed_umi >= thresh_a) %>% pull(cell_barcode)
  pos_i_cells <- crispri_df %>% filter(fixed_umi >= thresh_i) %>% pull(cell_barcode)

  message(paste("CRISPRa positive cells:", length(pos_a_cells)))
  message(paste("CRISPRi positive cells:", length(pos_i_cells)))

  # Variable guide assignment for positive cells
  # Rule: top guide >= 10 UMIs AND ratio (top / second) >= 5
  var_guides <- CIRI_long %>%
    filter(!feature %in% c(a_genes, i_genes)) %>%
    group_by(cell_barcode, type) %>%
    arrange(desc(umi), .by_group = TRUE) %>%
    summarise(
      top_guide   = first(feature),
      top_umi     = first(umi),
      second_umi  = ifelse(n() >= 2, nth(umi, 2), 0),
      ratio       = ifelse(second_umi > 0, top_umi / second_umi, Inf),
      assigned    = top_umi >= 10 & ratio >= 5,
      .groups     = "drop"
    )

  # Build annotation table
  var_a <- var_guides %>% filter(type == "a", assigned) %>%
    select(cell_barcode, feature_a = top_guide)
  var_i <- var_guides %>% filter(type == "i", assigned) %>%
    select(cell_barcode, feature_i = top_guide)

  annotation <- data.frame(cell_barcode = union(pos_a_cells, pos_i_cells)) %>%
    left_join(var_a, by = "cell_barcode") %>%
    left_join(var_i, by = "cell_barcode") %>%
    mutate(
      feature_a = ifelse(cell_barcode %in% pos_a_cells, replace_na(feature_a, "NA"), "NA"),
      feature_i = ifelse(cell_barcode %in% pos_i_cells, replace_na(feature_i, "NA"), "NA")
    )

  list(annotation = annotation, crispra = crispra_df, crispri = crispri_df)
}

# ---------------------------------------------------------------------------
# Assignment Logic - Dual Variable Guide (strategy = 2)
# ---------------------------------------------------------------------------

assign_dual <- function(CIRI_long, a_genes, i_genes, thresh_a, thresh_i) {
  # (Re-use CRISPRa/i fixed-guide summaries from single mode)
  base <- assign_single(CIRI_long, a_genes, i_genes, thresh_a, thresh_i)

  pos_a_cells <- base$annotation %>% filter(feature_a != "NA") %>% pull(cell_barcode)
  pos_i_cells <- base$annotation %>% filter(feature_i != "NA") %>% pull(cell_barcode)

  # Dual variable rule: sum of top 2 >= 4 UMIs AND ratio (top2_sum / third) >= 10
  # Cells with top two guides targeting DIFFERENT genes are removed.
  var_guides_dual <- CIRI_long %>%
    filter(!feature %in% c(a_genes, i_genes)) %>%
    group_by(cell_barcode, type) %>%
    arrange(desc(umi), .by_group = TRUE) %>%
    summarise(
      top_guide      = first(feature),
      second_guide   = ifelse(n() >= 2, nth(feature, 2), NA_character_),
      top_umi        = first(umi),
      second_umi     = ifelse(n() >= 2, nth(umi, 2), 0),
      third_umi      = ifelse(n() >= 3, nth(umi, 3), 0),
      top2_sum       = top_umi + second_umi,
      ratio          = ifelse(third_umi > 0, top2_sum / third_umi, Inf),
      same_gene      = !is.na(second_guide) &
                       (str_extract(top_guide, "^[^_]+") == str_extract(second_guide, "^[^_]+")),
      assigned       = top2_sum >= 4 & ratio >= 10 & same_gene,
      combined_guide = ifelse(assigned, paste(top_guide, second_guide, sep = ";"), NA_character_),
      .groups        = "drop"
    )

  var_a <- var_guides_dual %>% filter(type == "a", assigned) %>%
    select(cell_barcode, feature_a = combined_guide)
  var_i <- var_guides_dual %>% filter(type == "i", assigned) %>%
    select(cell_barcode, feature_i = combined_guide)

  annotation <- data.frame(cell_barcode = union(pos_a_cells, pos_i_cells)) %>%
    left_join(var_a, by = "cell_barcode") %>%
    left_join(var_i, by = "cell_barcode") %>%
    mutate(
      feature_a = replace_na(feature_a, "NA"),
      feature_i = replace_na(feature_i, "NA")
    )

  list(annotation = annotation, crispra = base$crispra, crispri = base$crispri)
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

main <- function() {
  params <- parse_args()

  message("=== CIRI Pipeline - Step 01: Perturbation Assignment ===")
  message(paste("Directory :", params$dir))
  message(paste("Matrix    :", params$matrix))
  message(paste("Strategy  :", params$strategy))

  # Load guide annotations
  ann_path <- file.path(params$dir, params$guides)
  if (!file.exists(ann_path)) stop(paste("guides.csv not found:", ann_path))

  ann <- read.csv(ann_path, header = FALSE)
  names(ann) <- c("feature", "type", "fixed")

  a_genes <- ann[ann$type == "a" & ann$fixed == "f", "feature"]
  i_genes <- ann[ann$type == "i" & ann$fixed == "f", "feature"]
  message(paste("Fixed CRISPRa guides:", paste(a_genes, collapse = ", ")))
  message(paste("Fixed CRISPRi guides:", paste(i_genes, collapse = ", ")))

  # Load CRISPR matrix
  h5_path <- file.path(params$dir, params$matrix)
  if (!file.exists(h5_path)) stop(paste("H5 matrix not found:", h5_path))

  message("Loading CRISPR capture data ...")
  crispr_mat <- load_crispr_h5(h5_path)

  CIRI_df   <- as.data.frame(t(as.matrix(crispr_mat)))
  CIRI_long <- CIRI_df %>%
    tibble::rownames_to_column("cell_barcode") %>%
    pivot_longer(cols = -cell_barcode, names_to = "feature", values_to = "umi") %>%
    left_join(ann, by = "feature") %>%
    group_by(cell_barcode) %>%
    mutate(
      total_umis = sum(umi),
      percentage = (umi / total_umis) * 100
    ) %>%
    ungroup()

  write.csv(CIRI_long, file.path(params$dir, "CIRI_long.csv"), row.names = FALSE)

  # Thresholds
  thresh_a <- if (params$threshold_a == -1) {
    umi_a <- CIRI_long %>% filter(type == "a", feature %in% a_genes) %>%
      group_by(cell_barcode) %>% summarise(s = sum(umi)) %>% pull(s)
    auto_threshold(umi_a, label = "CRISPRa")
  } else {
    params$threshold_a
  }

  thresh_i <- if (params$threshold_i == -1) {
    umi_i <- CIRI_long %>% filter(type == "i", feature %in% i_genes) %>%
      group_by(cell_barcode) %>% summarise(s = sum(umi)) %>% pull(s)
    auto_threshold(umi_i, label = "CRISPRi")
  } else {
    params$threshold_i
  }

  message(paste("Using threshold_a:", thresh_a))
  message(paste("Using threshold_i:", thresh_i))

  # Assign perturbations
  result <- if (params$strategy == 1L) {
    message("Running single variable guide assignment ...")
    assign_single(CIRI_long, a_genes, i_genes, thresh_a, thresh_i)
  } else {
    message("Running dual variable guide assignment ...")
    assign_dual(CIRI_long, a_genes, i_genes, thresh_a, thresh_i)
  }

  # Save QC plots
  message("Saving QC plots ...")
  save_qc_plots(CIRI_long, result$crispra, result$crispri, params$dir)

  # Save annotation
  out_ann <- file.path(params$dir, "annotation_data.csv")
  write.csv(result$annotation, out_ann, row.names = FALSE)
  message(paste("Annotation saved to:", out_ann))
  message(paste("Total assigned cells:", nrow(result$annotation)))

  message("Step 01 complete.")
}

main()
