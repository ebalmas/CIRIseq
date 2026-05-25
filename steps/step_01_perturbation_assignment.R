#!/usr/bin/env Rscript
# =============================================================================
# CIRI Pipeline — Step 01: Perturbation Assignment
# =============================================================================
# Reads the 10x H5 matrix, assigns CRISPRa / CRISPRi guide identities to
# each cell, and produces QC plots.
#
# INPUT  (provide directly):
#   --data_dir   folder containing the H5 matrix and guides.csv
#
# OUTPUT STRUCTURE:
#   Output/
#     <YYMMDD>_step01_assignment_<sample>/
#       csv/        CIRI_long.csv
#       plots/      total_umixguide.pdf, UMI histograms, fixed_vs_variable.pdf
#                   single_guide_plots/
#       stats/      assignment_summary.txt
#       R_objects/  (none — data passed as CSV)
#       to_scratch/ annotation_data.csv   ← consumed by Step 02
#
# Usage:
#   Rscript steps/step_01_perturbation_assignment.R \
#     --data_dir    /path/to/data \
#     --matrix      filtered_feature_bc_matrix.h5 \
#     --sample      AB011 \
#     --strategy    1
#   Rscript steps/step_01_perturbation_assignment.R --help
# =============================================================================

STEPS_DIR <- dirname(sys.frame(1)$ofile)
source(file.path(STEPS_DIR, "utils.R"))

SPEC <- list(
  data_dir    = list(flag="--data_dir",    required=TRUE,  type="character",
                     help="Directory containing the H5 matrix and guides.csv"),
  matrix      = list(flag="--matrix",      required=TRUE,  type="character",
                     help="H5 matrix filename inside --data_dir"),
  guides      = list(flag="--guides",      required=FALSE, type="character",
                     default="guides.csv",
                     help="Guide annotation CSV inside --data_dir"),
  output_root = list(flag="--output_root", required=FALSE, type="character",
                     default="Output",
                     help="Top-level output directory"),
  sample      = list(flag="--sample",      required=FALSE, type="character",
                     default="CIRI",
                     help="Sample / experiment name"),
  strategy    = list(flag="--strategy",    required=FALSE, type="integer",
                     default=1L,
                     help="1 = single variable guide, 2 = dual variable guides"),
  threshold_a = list(flag="--threshold_a", required=FALSE, type="numeric",
                     default=-1,
                     help="CRISPRa UMI threshold; -1 = auto-detect via KDE"),
  threshold_i = list(flag="--threshold_i", required=FALSE, type="numeric",
                     default=-1,
                     help="CRISPRi UMI threshold; -1 = auto-detect via KDE")
)

# ---- helpers ---------------------------------------------------------------
load_libs <- function() suppressMessages({
  library(dplyr); library(tidyr); library(stringr); library(ggplot2)
  library(quantmod); library(pracma); library(hdf5r)
  library(Matrix); library(zoo); library(scales)
})

auto_threshold <- function(umi_vec, label = "") {
  umi_vec <- umi_vec[umi_vec > 0]
  if (length(umi_vec) < 20) { log_warn("Too few UMIs for auto-threshold (", label, "). Using 5."); return(5) }
  dens   <- density(umi_vec, bw = "SJ", n = 1024)
  dy     <- diff(dens$y)
  valley <- which(dy[-length(dy)] < 0 & dy[-1] >= 0)
  if (length(valley) == 0) { log_warn("No KDE valley found (", label, "). Using median/2."); return(median(umi_vec)/2) }
  thresh <- dens$x[valley[1] + 1]
  log_info("Auto-threshold [", label, "]: ", round(thresh, 2))
  thresh
}

load_crispr_h5 <- function(h5_path) {
  f    <- H5File$new(h5_path, mode = "r"); on.exit(f$close_all(), add = TRUE)
  feat <- f[["matrix/features"]]
  idx  <- which(feat[["feature_type"]][] == "CRISPR Guide Capture")
  if (!length(idx)) stop("No CRISPR Guide Capture features found.")
  mat  <- new("dgCMatrix", x = as.numeric(f[["matrix/data"]][]),
              i = f[["matrix/indices"]][], p = f[["matrix/indptr"]][],
              Dim = f[["matrix/shape"]][])
  rownames(mat) <- feat[["name"]][]; colnames(mat) <- f[["matrix/barcodes"]][]
  mat[idx, , drop = FALSE]
}

fixed_summaries <- function(CIRI_long, a_genes, i_genes, thresh_a, thresh_i) {
  make_summary <- function(type_filter, fixed_genes, label) {
    CIRI_long %>% filter(type == type_filter) %>%
      mutate(is_fixed = feature %in% fixed_genes) %>%
      group_by(cell_barcode) %>%
      summarise(total_umi = sum(umi), fixed_umi = sum(umi[is_fixed]),
                percentage = (fixed_umi / total_umi) * 100,
                feature = label, .groups = "drop") %>%
      filter(total_umi > 0)
  }
  crispra <- make_summary("a", a_genes, "CRISPRa")
  crispri <- make_summary("i", i_genes, "CRISPRi")
  pos_a   <- filter(crispra, fixed_umi >= thresh_a) %>% pull(cell_barcode)
  pos_i   <- filter(crispri, fixed_umi >= thresh_i) %>% pull(cell_barcode)
  log_info("CRISPRa positive: ", length(pos_a), " cells")
  log_info("CRISPRi positive: ", length(pos_i), " cells")
  list(crispra = crispra, crispri = crispri, pos_a = pos_a, pos_i = pos_i)
}

assign_single <- function(CIRI_long, a_genes, i_genes, thresh_a, thresh_i) {
  fs  <- fixed_summaries(CIRI_long, a_genes, i_genes, thresh_a, thresh_i)
  var <- CIRI_long %>% filter(!feature %in% c(a_genes, i_genes)) %>%
    group_by(cell_barcode, type) %>% arrange(desc(umi), .by_group = TRUE) %>%
    summarise(top_guide = first(feature), top_umi = first(umi),
              second_umi = ifelse(n() >= 2, nth(umi, 2), 0),
              ratio = ifelse(second_umi > 0, top_umi / second_umi, Inf),
              assigned = top_umi >= 10 & ratio >= 5, .groups = "drop")
  var_a <- filter(var, type == "a", assigned) %>% select(cell_barcode, feature_a = top_guide)
  var_i <- filter(var, type == "i", assigned) %>% select(cell_barcode, feature_i = top_guide)
  ann <- data.frame(cell_barcode = union(fs$pos_a, fs$pos_i)) %>%
    left_join(var_a, by = "cell_barcode") %>%
    left_join(var_i, by = "cell_barcode") %>%
    mutate(feature_a = ifelse(cell_barcode %in% fs$pos_a, replace_na(feature_a, "NA"), "NA"),
           feature_i = ifelse(cell_barcode %in% fs$pos_i, replace_na(feature_i, "NA"), "NA"))
  list(annotation = ann, crispra = fs$crispra, crispri = fs$crispri)
}

assign_dual <- function(CIRI_long, a_genes, i_genes, thresh_a, thresh_i) {
  fs  <- fixed_summaries(CIRI_long, a_genes, i_genes, thresh_a, thresh_i)
  var <- CIRI_long %>% filter(!feature %in% c(a_genes, i_genes)) %>%
    group_by(cell_barcode, type) %>% arrange(desc(umi), .by_group = TRUE) %>%
    summarise(top_guide = first(feature), second_guide = ifelse(n()>=2, nth(feature,2), NA_character_),
              top_umi = first(umi), second_umi = ifelse(n()>=2, nth(umi,2), 0),
              third_umi = ifelse(n()>=3, nth(umi,3), 0),
              top2_sum = top_umi + second_umi,
              ratio = ifelse(third_umi > 0, top2_sum / third_umi, Inf),
              same_gene = !is.na(second_guide) &
                str_extract(top_guide,"^[^_]+") == str_extract(second_guide,"^[^_]+"),
              assigned = top2_sum >= 4 & ratio >= 10 & same_gene,
              combined = ifelse(assigned, paste(top_guide, second_guide, sep=";"), NA_character_),
              .groups = "drop")
  var_a <- filter(var, type=="a", assigned) %>% select(cell_barcode, feature_a = combined)
  var_i <- filter(var, type=="i", assigned) %>% select(cell_barcode, feature_i = combined)
  ann <- data.frame(cell_barcode = union(fs$pos_a, fs$pos_i)) %>%
    left_join(var_a, by="cell_barcode") %>% left_join(var_i, by="cell_barcode") %>%
    mutate(feature_a = replace_na(feature_a, "NA"), feature_i = replace_na(feature_i, "NA"))
  list(annotation = ann, crispra = fs$crispra, crispri = fs$crispri)
}

save_qc_plots <- function(CIRI_long, crispra, crispri, plots_dir) {
  plot_df <- CIRI_long %>%
    group_by(feature, type) %>% summarise(total_umi = sum(umi), .groups="drop") %>%
    arrange(type, feature) %>% mutate(feature = factor(feature, levels = unique(feature)))
  ggsave(file.path(plots_dir, "total_umixguide.pdf"),
    ggplot(plot_df, aes(x=feature, y=total_umi, fill=type)) +
      geom_bar(stat="identity") + theme_minimal() +
      labs(title="Total UMIs per guide") +
      theme(axis.text.x = element_text(angle=90, hjust=1)))

  umis_cell <- CIRI_long %>% group_by(cell_barcode) %>% summarise(total=sum(umi), .groups="drop")
  ggsave(file.path(plots_dir, "UMIxcell_hist.pdf"),
    ggplot(umis_cell, aes(x=total)) + geom_histogram(binwidth=1) + labs(title="Total UMIs per cell"))

  CRISPRai_plot <- bind_rows(crispra, crispri) %>% mutate(variable_umi = total_umi - fixed_umi)
  ggsave(file.path(plots_dir, "fixed_vs_variable_scatter.pdf"),
    ggplot(CRISPRai_plot, aes(x=variable_umi, y=fixed_umi, color=feature)) +
      geom_point(alpha=0.1, size=0.5) + theme_minimal() +
      labs(title="Fixed vs Variable UMI"),
    width=12, height=10)

  sg_dir <- file.path(plots_dir, "single_guide_plots")
  dir.create(sg_dir, showWarnings=FALSE)
  for (guide in unique(CIRI_long$feature)) {
    df <- dplyr::filter(CIRI_long, feature == guide)
    ggsave(file.path(sg_dir, paste0("umi_", guide, ".pdf")),
           ggplot(df, aes(x=umi)) + geom_histogram(binwidth=1) + labs(title=guide))
  }
}

# ---- main ------------------------------------------------------------------
main <- function() {
  load_libs()
  p   <- parse_args(SPEC)
  out <- make_out_dirs(p$output_root, "step01_assignment", p$sample)

  step_banner("01", "Perturbation Assignment",
    inputs  = c(file.path(p$data_dir, p$matrix), file.path(p$data_dir, p$guides)),
    outputs = c(out$csv, out$plots, out$stats,
                file.path(out$to_scratch, "annotation_data.csv"))
  )

  p$data_dir <- normalizePath(p$data_dir, mustWork = TRUE)
  assert_file(file.path(p$data_dir, p$matrix))
  assert_file(file.path(p$data_dir, p$guides))
  if (!p$strategy %in% c(1L, 2L)) stop("--strategy must be 1 or 2")

  ann <- read.csv(file.path(p$data_dir, p$guides), header=FALSE)
  names(ann) <- c("feature","type","fixed")
  a_genes <- ann[ann$type=="a" & ann$fixed=="f", "feature"]
  i_genes <- ann[ann$type=="i" & ann$fixed=="f", "feature"]
  log_info("Fixed CRISPRa guides: ", paste(a_genes, collapse=", "))
  log_info("Fixed CRISPRi guides: ", paste(i_genes, collapse=", "))

  log_info("Loading CRISPR capture data ...")
  crispr_mat <- load_crispr_h5(file.path(p$data_dir, p$matrix))
  CIRI_long  <- as.data.frame(t(as.matrix(crispr_mat))) %>%
    tibble::rownames_to_column("cell_barcode") %>%
    pivot_longer(-cell_barcode, names_to="feature", values_to="umi") %>%
    left_join(ann, by="feature") %>%
    group_by(cell_barcode) %>%
    mutate(total_umis=sum(umi), percentage=(umi/total_umis)*100) %>% ungroup()

  # csv/
  write.csv(CIRI_long, file.path(out$csv, "CIRI_long.csv"), row.names=FALSE)

  # thresholds
  thresh_a <- if (p$threshold_a == -1) {
    umi_a <- CIRI_long %>% filter(type=="a", feature %in% a_genes) %>%
      group_by(cell_barcode) %>% summarise(s=sum(umi), .groups="drop") %>% pull(s)
    auto_threshold(umi_a, "CRISPRa")
  } else p$threshold_a

  thresh_i <- if (p$threshold_i == -1) {
    umi_i <- CIRI_long %>% filter(type=="i", feature %in% i_genes) %>%
      group_by(cell_barcode) %>% summarise(s=sum(umi), .groups="drop") %>% pull(s)
    auto_threshold(umi_i, "CRISPRi")
  } else p$threshold_i

  result <- if (p$strategy == 1L) assign_single(CIRI_long, a_genes, i_genes, thresh_a, thresh_i) \
            else                  assign_dual(CIRI_long,   a_genes, i_genes, thresh_a, thresh_i)

  # plots/
  log_info("Saving QC plots ...")
  save_qc_plots(CIRI_long, result$crispra, result$crispri, out$plots)

  # stats/
  n_a  <- sum(result$annotation$feature_a != "NA")
  n_i  <- sum(result$annotation$feature_i != "NA")
  stats_txt <- c(
    paste("Sample          :", p$sample),
    paste("Strategy        :", p$strategy),
    paste("Threshold CRISPRa:", thresh_a),
    paste("Threshold CRISPRi:", thresh_i),
    paste("Total assigned cells :", nrow(result$annotation)),
    paste("  with CRISPRa guide :", n_a),
    paste("  with CRISPRi guide :", n_i)
  )
  writeLines(stats_txt, file.path(out$stats, "assignment_summary.txt"))
  message(paste(stats_txt, collapse="\n"))

  # to_scratch/ — annotation consumed by Step 02
  write.csv(result$annotation,
            file.path(out$to_scratch, "annotation_data.csv"), row.names=FALSE)

  log_info("Step 01 complete.")
  log_info("Output folder: ", out$out_dir)
}

main()
