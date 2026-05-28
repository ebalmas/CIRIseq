# =============================================================================
# CIRI — Step 01: Perturbation Assignment
# =============================================================================

#' Assign CRISPR guide identities to cells
#'
#' Reads either a 10x Genomics H5 matrix or the CellRanger
#' \code{protospacer_calls_per_cell.csv} file, matches guides against the
#' annotation CSV, and assigns CRISPRa / CRISPRi identities to each cell.
#'
#' @section Guide name matching:
#' CellRanger collapses guide replicates (e.g. \code{ATF7IP_1A} and
#' \code{ATF7IP_1B} both become \code{ATF7IP_1} in the H5 and the
#' protospacer CSV). The function automatically harmonises names from
#' \code{guides} against the names in the data using a three-step lookup:
#' strip trailing A/B → gene name (before first \code{_}) → exact match.
#' A warning is printed for any guide in \code{guides} that cannot be matched.
#'
#' @param data_dir    Directory containing the input files.
#' @param matrix      H5 matrix filename inside \code{data_dir}. Used only when
#'                    \code{protospacer} is \code{NULL}.
#' @param protospacer \code{protospacer_calls_per_cell.csv} filename inside
#'                    \code{data_dir}. When provided this is used instead of
#'                    the H5 (recommended — avoids HDF5 name mismatches).
#' @param guides      Guide annotation CSV filename inside \code{data_dir}.
#'                    Format (no header): \code{feature, type (a/i), fixed (f/v)}.
#'                    Default \code{"guides.csv"}.
#' @param sample      Experiment name, used in output folder name.
#' @param output_root Top-level output directory. Default \code{"Output"}.
#' @param strategy    Integer. \code{1} = single arm: cells assigned if they have
#'                    a valid CRISPRa OR CRISPRi variable guide (or both).
#'                    \code{2} = CIRI dual arm: cells assigned only if they have
#'                    a valid variable guide for BOTH CRISPRa AND CRISPRi.
#'                    Use \code{2} for CRISPRa+CRISPRi combinatorial screens.
#' @param threshold_a CRISPRa fixed-guide UMI threshold. \code{-1} = auto-detect (default).
#' @param threshold_i CRISPRi fixed-guide UMI threshold. \code{-1} = auto-detect (default).
#'
#' @return Invisibly returns the path to the output folder.
#'
#' @section Output structure:
#' \preformatted{
#' Output/
#'   <YYMMDD>_step01_assignment_<sample>/
#'     csv/        CIRI_long.csv, guide_name_mapping.csv
#'     plots/      total_umixguide.pdf, fixed_vs_variable_scatter.pdf,
#'                 UMIxcell_hist.pdf, single_guide_plots/
#'     stats/      assignment_summary.txt
#'     to_scratch/ annotation_data.csv
#' }
#'
#' @examples
#' \dontrun{
#' # Recommended: use protospacer_calls_per_cell.csv
#' ciri_step01_assignment(
#'   data_dir    = "/data/AB011",
#'   protospacer = "protospacer_calls_per_cell.csv",
#'   guides      = "guides_2.csv",
#'   sample      = "AB011",
#'   strategy    = 2
#' )
#'
#' # Alternative: use H5 matrix directly
#' ciri_step01_assignment(
#'   data_dir = "/data/AB011",
#'   matrix   = "filtered_feature_bc_matrix.h5",
#'   guides   = "guides_2.csv",
#'   sample   = "AB011",
#'   strategy = 2
#' )
#' }
#'
#' @importFrom dplyr filter mutate group_by summarise left_join arrange pull
#'   bind_rows select ungroup rename
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom ggplot2 ggplot aes geom_bar geom_histogram geom_point
#'   theme_minimal labs ggsave scale_y_continuous theme element_text
#' @importFrom scales label_number
#' @importFrom utils write.csv read.csv
#' @export
ciri_step01_assignment <- function(data_dir,
                                   matrix       = NULL,
                                   protospacer  = NULL,
                                   guides       = "guides.csv",
                                   sample       = "CIRI",
                                   output_root  = "Output",
                                   strategy     = 1L,
                                   threshold_a  = -1,
                                   threshold_i  = -1) {

  if (is.null(matrix) && is.null(protospacer))
    stop("Provide either 'matrix' (H5 filename) or 'protospacer' (protospacer_calls_per_cell.csv).",
         call. = FALSE)

  out      <- make_out_dirs(output_root, "step01_assignment", sample)
  data_dir <- normalizePath(data_dir, mustWork = TRUE)

  input_file <- if (!is.null(protospacer))
    file.path(data_dir, protospacer) else file.path(data_dir, matrix)

  step_banner("01", "Perturbation Assignment",
    inputs  = c(input_file, file.path(data_dir, guides)),
    outputs = c(out$csv, out$plots, out$stats,
                file.path(out$to_scratch, "annotation_data.csv"))
  )

  assert_file(input_file)
  assert_file(file.path(data_dir, guides))
  if (!strategy %in% c(1L, 2L)) stop("strategy must be 1 or 2", call. = FALSE)

  # ---- load and harmonise guide annotations --------------------------------
  ann_raw <- read.csv(file.path(data_dir, guides), header = FALSE)
  names(ann_raw) <- c("feature", "type", "fixed")

  # ---- load CRISPR data ----------------------------------------------------
  if (!is.null(protospacer)) {
    log_info("Reading protospacer_calls_per_cell.csv ...")
    result <- .load_from_protospacer(file.path(data_dir, protospacer), ann_raw)
  } else {
    .check_cran_pkgs("hdf5r")
    log_info("Loading CRISPR capture data from H5 ...")
    result <- .load_from_h5(file.path(data_dir, matrix), ann_raw)
  }

  CIRI_long <- result$CIRI_long
  ann       <- result$ann          # harmonised annotation (matched names)
  mapping   <- result$mapping      # name mapping table for reference

  # Save name mapping so user can inspect
  write.csv(mapping, file.path(out$csv, "guide_name_mapping.csv"), row.names = FALSE)
  log_info("Guide name mapping saved to csv/guide_name_mapping.csv")

  a_genes <- ann[ann$type == "a" & ann$fixed == "f", "feature"]
  i_genes <- ann[ann$type == "i" & ann$fixed == "f", "feature"]
  log_info("Fixed CRISPRa: ", paste(a_genes, collapse = ", "))
  log_info("Fixed CRISPRi: ", paste(i_genes, collapse = ", "))

  write.csv(CIRI_long, file.path(out$csv, "CIRI_long.csv"), row.names = FALSE)

  # ---- thresholds ----------------------------------------------------------
  # KDE diagnostic plots saved to plots/ so you can verify the cut visually
  thresh_a <- .get_threshold(CIRI_long, "a", a_genes, threshold_a, "CRISPRa",
                              plots_dir = out$plots)
  thresh_i <- .get_threshold(CIRI_long, "i", i_genes, threshold_i, "CRISPRi",
                              plots_dir = out$plots)
  log_info("Threshold CRISPRa: ", thresh_a, " | CRISPRi: ", thresh_i)

  # ---- assign --------------------------------------------------------------
  assign_result <- if (strategy == 1L)
    .assign_single(CIRI_long, a_genes, i_genes, thresh_a, thresh_i)
  else
    .assign_dual(CIRI_long, a_genes, i_genes, thresh_a, thresh_i)

  # ---- plots ---------------------------------------------------------------
  log_info("Saving QC plots ...")
  .save_assignment_plots(CIRI_long, assign_result$crispra, assign_result$crispri, out$plots)

  # ---- stats ---------------------------------------------------------------
  stats_txt <- c(
    paste("Sample            :", sample),
    paste("Strategy          :", strategy),
    paste("Input             :", if (!is.null(protospacer)) "protospacer_calls_per_cell.csv" else "H5 matrix"),
    paste("Threshold CRISPRa :", thresh_a),
    paste("Threshold CRISPRi :", thresh_i),
    paste("Assigned cells    :", nrow(assign_result$annotation)),
    paste("  with CRISPRa    :", sum(assign_result$annotation$feature_a != "NA")),
    paste("  with CRISPRi    :", sum(assign_result$annotation$feature_i != "NA"))
  )
  writeLines(stats_txt, file.path(out$stats, "assignment_summary.txt"))
  message(paste(stats_txt, collapse = "\n"))

  write.csv(assign_result$annotation,
            file.path(out$to_scratch, "annotation_data.csv"),
            row.names = FALSE, quote = TRUE)

  log_info("Step 01 complete. Output: ", out$out_dir)
  invisible(out$out_dir)
}

# =============================================================================
# Private helpers
# =============================================================================

# ---------------------------------------------------------------------------
# .harmonise_guide_names
# Maps guides.csv feature names to the (possibly collapsed) names used by
# CellRanger in the H5 or protospacer CSV.
# Returns: ann data.frame with 'feature' column replaced by the matched name,
#          plus a mapping table for reference.
# ---------------------------------------------------------------------------
.harmonise_guide_names <- function(ann_raw, data_guide_names) {
  find_match <- function(csv_name, data_names) {
    # 1. strip trailing A or B after _number: "ATF7IP_1A" -> "ATF7IP_1"
    m <- regmatches(csv_name, regexpr("^.+_\\d+(?=[AB]$)", csv_name, perl = TRUE))
    if (length(m) && m %in% data_names) return(m)
    # 2. gene name only (before first underscore): "NTCa_1A" -> "NTCa"
    gene <- strsplit(csv_name, "_")[[1]][1]
    if (gene %in% data_names) return(gene)
    # 3. exact match
    if (csv_name %in% data_names) return(csv_name)
    return(NA_character_)
  }

  ann_raw$matched <- vapply(ann_raw$feature, find_match,
                            character(1), data_names = data_guide_names)

  unmatched <- ann_raw$feature[is.na(ann_raw$matched)]
  if (length(unmatched))
    log_warn(length(unmatched), " guide(s) in guides CSV not found in data: ",
             paste(unmatched, collapse = ", "))

  # Deduplicate: multiple csv rows can map to same data name (A/B replicates)
  # Keep first occurrence for type/fixed classification
  ann_matched <- ann_raw[!is.na(ann_raw$matched), ]
  ann_dedup   <- ann_matched[!duplicated(ann_matched$matched), ]
  ann_out     <- data.frame(feature = ann_dedup$matched,
                            type    = ann_dedup$type,
                            fixed   = ann_dedup$fixed,
                            stringsAsFactors = FALSE)

  mapping <- data.frame(
    guides_csv_name    = ann_raw$feature,
    cellranger_name    = ann_raw$matched,
    type               = ann_raw$type,
    fixed              = ann_raw$fixed,
    stringsAsFactors   = FALSE
  )

  list(ann = ann_out, mapping = mapping)
}

# ---------------------------------------------------------------------------
# .load_from_protospacer
# Reads protospacer_calls_per_cell.csv and builds CIRI_long.
# ---------------------------------------------------------------------------
.load_from_protospacer <- function(proto_path, ann_raw) {
  proto <- read.csv(proto_path, header = TRUE, stringsAsFactors = FALSE)
  # Expected columns: cell_barcode, num_features, feature_call, num_umis

  # Parse pipe-separated guide calls into long format
  rows <- lapply(seq_len(nrow(proto)), function(i) {
    features <- strsplit(proto$feature_call[i], "\\|")[[1]]
    umis     <- as.integer(strsplit(as.character(proto$num_umis[i]), "\\|")[[1]])
    if (length(features) == 0 || all(is.na(features))) return(NULL)
    data.frame(cell_barcode = proto$cell_barcode[i],
               feature      = features,
               umi          = umis,
               stringsAsFactors = FALSE)
  })
  CIRI_raw <- do.call(rbind, rows[!sapply(rows, is.null)])

  # Harmonise guide names
  data_guide_names <- unique(CIRI_raw$feature)
  harm <- .harmonise_guide_names(ann_raw, data_guide_names)

  # Join annotation
  CIRI_long <- dplyr::left_join(CIRI_raw, harm$ann, by = "feature") |>
    dplyr::group_by(cell_barcode) |>
    dplyr::mutate(total_umis   = sum(umi),
                  percentage   = umi / total_umis * 100) |>
    dplyr::ungroup()

  list(CIRI_long = CIRI_long, ann = harm$ann, mapping = harm$mapping)
}

# ---------------------------------------------------------------------------
# .load_from_h5
# Reads CRISPR capture rows from H5 and builds CIRI_long.
# ---------------------------------------------------------------------------
.load_from_h5 <- function(h5_path, ann_raw) {
  crispr_mat <- .load_crispr_h5(h5_path)

  CIRI_raw <- as.data.frame(t(as.matrix(crispr_mat))) |>
    tibble::rownames_to_column("cell_barcode") |>
    tidyr::pivot_longer(-cell_barcode, names_to = "feature", values_to = "umi")

  # Harmonise guide names
  data_guide_names <- unique(CIRI_raw$feature)
  harm <- .harmonise_guide_names(ann_raw, data_guide_names)

  CIRI_long <- dplyr::left_join(CIRI_raw, harm$ann, by = "feature") |>
    dplyr::group_by(cell_barcode) |>
    dplyr::mutate(total_umis = sum(umi),
                  percentage = umi / total_umis * 100) |>
    dplyr::ungroup()

  list(CIRI_long = CIRI_long, ann = harm$ann, mapping = harm$mapping)
}

.load_crispr_h5 <- function(h5_path) {
  f    <- hdf5r::H5File$new(h5_path, mode = "r")
  on.exit(f$close_all(), add = TRUE)
  feat <- f[["matrix/features"]]
  idx  <- which(feat[["feature_type"]][] == "CRISPR Guide Capture")
  if (!length(idx)) stop("No CRISPR Guide Capture features found in H5.", call. = FALSE)
  mat  <- methods::new("dgCMatrix",
                       x   = as.numeric(f[["matrix/data"]][]),
                       i   = f[["matrix/indices"]][],
                       p   = f[["matrix/indptr"]][],
                       Dim = f[["matrix/shape"]][])
  rownames(mat) <- feat[["name"]][]
  colnames(mat) <- f[["matrix/barcodes"]][]
  mat[idx, , drop = FALSE]
}

# ---------------------------------------------------------------------------
# .auto_threshold
# Finds the threshold between "cells that did NOT receive the fixed guide"
# (noise peak, low UMIs) and "cells that DID receive it" (signal peak, high
# UMIs) using the first valley in the KDE of per-cell fixed-guide UMI sums.
#
# What the plot shows (saved to plots/threshold_kde_<label>.pdf):
#   • Grey histogram  — raw UMI distribution (log10 x-axis)
#   • Blue curve      — kernel density estimate (KDE)
#   • All valleys     — grey dotted vertical lines (all local minima found)
#   • Red dashed line — the chosen threshold (first/leftmost valley)
#   • Annotation      — threshold value printed on the plot
#
# If the red line lands in the wrong place, pass the threshold manually:
#   ciri_step01_assignment(..., threshold_a = 50, threshold_i = 30)
# ---------------------------------------------------------------------------
.auto_threshold <- function(umi_vec, label = "", plots_dir = NULL) {
  umi_vec <- umi_vec[umi_vec > 0]
  if (length(umi_vec) < 20) {
    log_warn("Too few UMIs for auto-threshold (", label, "). Using 5.")
    return(5)
  }

  dens   <- stats::density(umi_vec, bw = "SJ", n = 1024)
  dy     <- diff(dens$y)

  # All local minima (valleys): derivative goes from negative to positive
  all_valleys <- which(dy[-length(dy)] < 0 & dy[-1] >= 0)

  if (!length(all_valleys)) {
    log_warn("No KDE valley found (", label, "). Using median/2. ",
             "Check threshold_kde_", label, ".pdf and set threshold manually.")
    thresh <- stats::median(umi_vec) / 2
  } else {
    thresh <- dens$x[all_valleys[1] + 1]
  }

  log_info("Auto-threshold [", label, "]: ", round(thresh, 2),
           "  (", sum(umi_vec >= thresh), " cells above threshold)")

  # --- Save diagnostic KDE plot -------------------------------------------
  if (!is.null(plots_dir)) {
    df_dens <- data.frame(x = dens$x, y = dens$y)

    p <- ggplot2::ggplot() +
      # Histogram of UMIs (log10 scale)
      ggplot2::geom_histogram(
        data    = data.frame(umi = umi_vec),
        mapping = ggplot2::aes(x = umi, y = ggplot2::after_stat(density)),
        bins    = 80, fill = "grey80", colour = "grey60", alpha = 0.6
      ) +
      # KDE curve
      ggplot2::geom_line(
        data    = df_dens,
        mapping = ggplot2::aes(x = x, y = y),
        colour  = "#2166AC", linewidth = 1
      ) +
      # All valleys (grey dotted)
      {if (length(all_valleys))
        ggplot2::geom_vline(
          xintercept = dens$x[all_valleys + 1],
          linetype = "dotted", colour = "grey40", linewidth = 0.7
        )
      } +
      # Chosen threshold (red dashed)
      ggplot2::geom_vline(
        xintercept = thresh,
        linetype = "dashed", colour = "#D6604D", linewidth = 1.2
      ) +
      # Threshold label
      ggplot2::annotate(
        "text", x = thresh, y = Inf, vjust = 2, hjust = -0.15,
        label = paste0("threshold = ", round(thresh, 1)),
        colour = "#D6604D", size = 4
      ) +
      ggplot2::scale_x_log10() +
      ggplot2::labs(
        title    = paste0("Fixed-guide UMI threshold — ", label),
        subtitle = paste0(
          "KDE valley method  |  threshold = ", round(thresh, 1),
          "  |  cells above = ", sum(umi_vec >= thresh),
          " / ", length(umi_vec),
          "\nIf the red line looks wrong, set threshold_",
          tolower(substr(label, nchar(label), nchar(label))),
          " manually in ciri_step01_assignment()"
        ),
        x = "Fixed-guide UMIs per cell  (log10)",
        y = "Density"
      ) +
      ggplot2::theme_bw(base_size = 13) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.subtitle    = ggplot2::element_text(size = 9, colour = "grey40")
      )

    out_pdf <- file.path(plots_dir,
                         paste0("threshold_kde_", gsub(" ","_",label), ".pdf"))
    ggplot2::ggsave(out_pdf, p, width = 7, height = 5)
    log_info("KDE threshold plot saved: plots/", basename(out_pdf))
  }

  thresh
}

# Wrapper that pulls the right UMI vector and calls .auto_threshold
.get_threshold <- function(CIRI_long, type_filter, fixed_genes,
                            threshold, label, plots_dir = NULL) {
  if (threshold != -1) {
    log_info("Manual threshold [", label, "]: ", threshold)
    return(threshold)
  }
  umi_vec <- CIRI_long |>
    dplyr::filter(type == type_filter, feature %in% fixed_genes) |>
    dplyr::group_by(cell_barcode) |>
    dplyr::summarise(s = sum(umi), .groups = "drop") |>
    dplyr::pull(s)
  .auto_threshold(umi_vec, label, plots_dir = plots_dir)
}

.fixed_summaries <- function(CIRI_long, a_genes, i_genes, thresh_a, thresh_i) {
  make_s <- function(type_f, fixed_g, label) {
    CIRI_long |> dplyr::filter(type == type_f) |>
      dplyr::mutate(is_fixed = feature %in% fixed_g) |>
      dplyr::group_by(cell_barcode) |>
      dplyr::summarise(total_umi  = sum(umi),
                       fixed_umi  = sum(umi[is_fixed]),
                       percentage = fixed_umi / total_umi * 100,
                       feature    = label,
                       .groups    = "drop") |>
      dplyr::filter(total_umi > 0)
  }
  crispra <- make_s("a", a_genes, "CRISPRa")
  crispri <- make_s("i", i_genes, "CRISPRi")
  pos_a   <- dplyr::filter(crispra, fixed_umi >= thresh_a) |> dplyr::pull(cell_barcode)
  pos_i   <- dplyr::filter(crispri, fixed_umi >= thresh_i) |> dplyr::pull(cell_barcode)
  log_info("CRISPRa positive: ", length(pos_a), " | CRISPRi positive: ", length(pos_i))
  list(crispra = crispra, crispri = crispri, pos_a = pos_a, pos_i = pos_i)
}

# .assign_var_guides: assigns one variable guide per cell for one type arm.
# A cell is assigned if: top guide UMI >= min_umi AND top/second ratio >= min_ratio.
.assign_var_guides <- function(CIRI_long, fixed_genes, type_filter,
                               min_umi = 10, min_ratio = 5) {
  CIRI_long |>
    dplyr::filter(type == type_filter, !feature %in% fixed_genes) |>
    dplyr::group_by(cell_barcode) |>
    dplyr::arrange(dplyr::desc(umi), .by_group = TRUE) |>
    dplyr::summarise(
      top_guide  = dplyr::first(feature),
      top_umi    = dplyr::first(umi),
      second_umi = ifelse(dplyr::n() >= 2, dplyr::nth(umi, 2), 0),
      ratio      = ifelse(second_umi > 0, top_umi / second_umi, Inf),
      assigned   = top_umi >= min_umi & ratio >= min_ratio,
      .groups    = "drop"
    ) |>
    dplyr::filter(assigned) |>
    dplyr::select(cell_barcode, top_guide)
}

# Strategy 1: cells assigned if they have a valid CRISPRa OR CRISPRi variable guide.
.assign_single <- function(CIRI_long, a_genes, i_genes, thresh_a, thresh_i) {
  fs    <- .fixed_summaries(CIRI_long, a_genes, i_genes, thresh_a, thresh_i)
  var_a <- .assign_var_guides(CIRI_long, c(a_genes, i_genes), "a") |>
    dplyr::rename(feature_a = top_guide)
  var_i <- .assign_var_guides(CIRI_long, c(a_genes, i_genes), "i") |>
    dplyr::rename(feature_i = top_guide)
  ann <- data.frame(cell_barcode = union(fs$pos_a, fs$pos_i)) |>
    dplyr::left_join(var_a, by = "cell_barcode") |>
    dplyr::left_join(var_i, by = "cell_barcode") |>
    dplyr::mutate(
      feature_a = ifelse(cell_barcode %in% fs$pos_a, tidyr::replace_na(feature_a, "NA"), "NA"),
      feature_i = ifelse(cell_barcode %in% fs$pos_i, tidyr::replace_na(feature_i, "NA"), "NA")
    )
  list(annotation = ann, crispra = fs$crispra, crispri = fs$crispri)
}

# Strategy 2 — CIRI dual arm: cells must have BOTH a CRISPRa AND CRISPRi variable guide.
.assign_dual <- function(CIRI_long, a_genes, i_genes, thresh_a, thresh_i) {
  fs    <- .fixed_summaries(CIRI_long, a_genes, i_genes, thresh_a, thresh_i)
  var_a <- .assign_var_guides(CIRI_long, c(a_genes, i_genes), "a") |>
    dplyr::rename(feature_a = top_guide)
  var_i <- .assign_var_guides(CIRI_long, c(a_genes, i_genes), "i") |>
    dplyr::rename(feature_i = top_guide)

  ciri_cells <- intersect(
    intersect(fs$pos_a, var_a$cell_barcode),
    intersect(fs$pos_i, var_i$cell_barcode)
  )
  only_a <- setdiff(intersect(fs$pos_a, var_a$cell_barcode), fs$pos_i)
  only_i <- setdiff(intersect(fs$pos_i, var_i$cell_barcode), fs$pos_a)
  log_info("CIRI (dual CRISPRa+CRISPRi) cells: ", length(ciri_cells))
  log_info("CRISPRa-only cells: ", length(only_a),
           " | CRISPRi-only cells: ", length(only_i))

  ann <- data.frame(cell_barcode = union(union(ciri_cells, only_a), only_i)) |>
    dplyr::left_join(var_a, by = "cell_barcode") |>
    dplyr::left_join(var_i, by = "cell_barcode") |>
    dplyr::mutate(
      feature_a = tidyr::replace_na(feature_a, "NA"),
      feature_i = tidyr::replace_na(feature_i, "NA")
    )
  list(annotation = ann, crispra = fs$crispra, crispri = fs$crispri)
}

.save_assignment_plots <- function(CIRI_long, crispra, crispri, plots_dir) {
  plot_df <- CIRI_long |>
    dplyr::filter(!is.na(type)) |>
    dplyr::group_by(feature, type) |>
    dplyr::summarise(total_umi = sum(umi), .groups = "drop") |>
    dplyr::arrange(type, feature) |>
    dplyr::mutate(feature = factor(feature, levels = unique(feature)))

  ggplot2::ggsave(file.path(plots_dir, "total_umixguide.pdf"),
    ggplot2::ggplot(plot_df, ggplot2::aes(x = feature, y = total_umi, fill = type)) +
      ggplot2::geom_bar(stat = "identity") + ggplot2::theme_minimal() +
      ggplot2::labs(title = "Total UMIs per guide") +
      ggplot2::scale_y_continuous(labels = scales::label_number()) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)))

  umis_cell <- CIRI_long |> dplyr::group_by(cell_barcode) |>
    dplyr::summarise(total = sum(umi), .groups = "drop")
  ggplot2::ggsave(file.path(plots_dir, "UMIxcell_hist.pdf"),
    ggplot2::ggplot(umis_cell, ggplot2::aes(x = total)) +
      ggplot2::geom_histogram(binwidth = 1) +
      ggplot2::labs(title = "UMIs per cell"))

  CRISPRai <- dplyr::bind_rows(crispra, crispri) |>
    dplyr::mutate(variable_umi = total_umi - fixed_umi)
  ggplot2::ggsave(file.path(plots_dir, "fixed_vs_variable_scatter.pdf"),
    ggplot2::ggplot(CRISPRai, ggplot2::aes(x = variable_umi, y = fixed_umi,
                                            color = feature)) +
      ggplot2::geom_point(alpha = 0.1, size = 0.5) + ggplot2::theme_minimal() +
      ggplot2::labs(title = "Fixed vs Variable UMI"),
    width = 12, height = 10)

  sg_dir <- file.path(plots_dir, "single_guide_plots")
  dir.create(sg_dir, showWarnings = FALSE)
  for (guide in unique(stats::na.omit(CIRI_long$feature))) {
    df <- dplyr::filter(CIRI_long, feature == guide)
    ggplot2::ggsave(file.path(sg_dir, paste0("umi_", guide, ".pdf")),
      ggplot2::ggplot(df, ggplot2::aes(x = umi)) +
        ggplot2::geom_histogram(binwidth = 1) +
        ggplot2::labs(title = guide))
  }
}

#' Assign CRISPR guide identities to cells
#'
#' Reads the 10x Genomics H5 matrix, identifies CRISPR Guide Capture features,
#' assigns CRISPRa / CRISPRi guide identities to each cell, and generates QC
#' plots. Thresholds can be set manually or auto-detected via KDE valley.
#'
#' @param data_dir    Directory containing the H5 matrix and guides CSV.
#' @param matrix      H5 matrix filename inside \code{data_dir}.
#' @param sample      Experiment name, used in output folder name.
#' @param output_root Top-level output directory. Default \code{"Output"}.
#' @param guides      Guide annotation CSV filename inside \code{data_dir}.
#'                    Format (no header): \code{feature, type (a/i), fixed (f/v)}.
#'                    Default \code{"guides.csv"}.
#' @param strategy    Integer. \code{1} = single arm: cells assigned if they have
#'                    a valid CRISPRa OR CRISPRi variable guide (or both).
#'                    \code{2} = CIRI dual arm: cells assigned only if they have
#'                    a valid variable guide for BOTH CRISPRa AND CRISPRi.
#'                    Use \code{2} for CRISPRa+CRISPRi combinatorial screens.
#' @param threshold_a CRISPRa UMI threshold. \code{-1} = auto-detect (default).
#' @param threshold_i CRISPRi UMI threshold. \code{-1} = auto-detect (default).
#'
#' @return Invisibly returns the path to the output folder.
#'
#' @section Output structure:
#' \preformatted{
#' Output/
#'   <YYMMDD>_step01_assignment_<sample>/
#'     csv/        CIRI_long.csv
#'     plots/      total_umixguide.pdf, fixed_vs_variable_scatter.pdf,
#'                 UMIxcell_hist.pdf, single_guide_plots/
#'     stats/      assignment_summary.txt
#'     to_scratch/ annotation_data.csv    <- read by Step 02 from scratch/
#' }
#'
#' @examples
#' \dontrun{
#' ciri_step01_assignment(
#'   data_dir = "/data/AB011",
#'   matrix   = "filtered_feature_bc_matrix.h5",
#'   sample   = "AB011",
#'   strategy = 1
#' )
#' }
#'
#' @importFrom dplyr filter mutate group_by summarise left_join arrange pull bind_rows
#' @importFrom tidyr pivot_longer replace_na
#' @importFrom stringr str_extract
#' @importFrom ggplot2 ggplot aes geom_bar geom_histogram geom_point theme_minimal
#'   labs ggsave scale_y_continuous theme element_text
#' @importFrom scales label_number
#' @importFrom utils write.csv read.csv
#' @export
ciri_step01_assignment <- function(data_dir,
                                   matrix,
                                   sample       = "CIRI",
                                   output_root  = "Output",
                                   guides       = "guides.csv",
                                   strategy     = 1L,
                                   threshold_a  = -1,
                                   threshold_i  = -1) {

  .check_cran_pkgs("hdf5r")

  out <- make_out_dirs(output_root, "step01_assignment", sample)
  data_dir <- normalizePath(data_dir, mustWork = TRUE)

  step_banner("01", "Perturbation Assignment",
    inputs  = c(file.path(data_dir, matrix), file.path(data_dir, guides)),
    outputs = c(out$csv, out$plots, out$stats,
                file.path(out$to_scratch, "annotation_data.csv"))
  )

  assert_file(file.path(data_dir, matrix))
  assert_file(file.path(data_dir, guides))
  if (!strategy %in% c(1L, 2L)) stop("strategy must be 1 or 2")

  # ---- guide annotations --------------------------------------------------
  ann <- read.csv(file.path(data_dir, guides), header = FALSE)
  names(ann) <- c("feature", "type", "fixed")
  a_genes <- ann[ann$type == "a" & ann$fixed == "f", "feature"]
  i_genes <- ann[ann$type == "i" & ann$fixed == "f", "feature"]
  log_info("Fixed CRISPRa: ", paste(a_genes, collapse = ", "))
  log_info("Fixed CRISPRi: ", paste(i_genes, collapse = ", "))

  # ---- load CRISPR matrix -------------------------------------------------
  log_info("Loading CRISPR capture data ...")
  crispr_mat <- .load_crispr_h5(file.path(data_dir, matrix))
  CIRI_long  <- as.data.frame(t(as.matrix(crispr_mat))) |>
    tibble::rownames_to_column("cell_barcode") |>
    tidyr::pivot_longer(-cell_barcode, names_to = "feature", values_to = "umi") |>
    dplyr::left_join(ann, by = "feature") |>
    dplyr::group_by(cell_barcode) |>
    dplyr::mutate(total_umis = sum(umi), percentage = umi / total_umis * 100) |>
    dplyr::ungroup()

  write.csv(CIRI_long, file.path(out$csv, "CIRI_long.csv"), row.names = FALSE)

  # ---- thresholds ---------------------------------------------------------
  thresh_a <- .get_threshold(CIRI_long, "a", a_genes, threshold_a, "CRISPRa",
                              plots_dir = out$plots)
  thresh_i <- .get_threshold(CIRI_long, "i", i_genes, threshold_i, "CRISPRi",
                              plots_dir = out$plots)
  log_info("Threshold CRISPRa: ", thresh_a, " | CRISPRi: ", thresh_i)

  # ---- assign -------------------------------------------------------------
  result <- if (strategy == 1L)
    .assign_single(CIRI_long, a_genes, i_genes, thresh_a, thresh_i)
  else
    .assign_dual(CIRI_long, a_genes, i_genes, thresh_a, thresh_i)

  # ---- plots --------------------------------------------------------------
  log_info("Saving QC plots ...")
  .save_assignment_plots(CIRI_long, result$crispra, result$crispri, out$plots)

  # ---- stats --------------------------------------------------------------
  stats_txt <- c(
    paste("Sample            :", sample),
    paste("Strategy          :", strategy),
    paste("Threshold CRISPRa :", thresh_a),
    paste("Threshold CRISPRi :", thresh_i),
    paste("Assigned cells    :", nrow(result$annotation)),
    paste("  with CRISPRa    :", sum(result$annotation$feature_a != "NA")),
    paste("  with CRISPRi    :", sum(result$annotation$feature_i != "NA"))
  )
  writeLines(stats_txt, file.path(out$stats, "assignment_summary.txt"))
  message(paste(stats_txt, collapse = "\n"))

  # ---- to_scratch/ --------------------------------------------------------
  write.csv(result$annotation,
            file.path(out$to_scratch, "annotation_data.csv"),
            row.names = FALSE, quote = TRUE)

  log_info("Step 01 complete. Output: ", out$out_dir)
  invisible(out$out_dir)
}

# ---- private helpers -------------------------------------------------------

.load_crispr_h5 <- function(h5_path) {
  f    <- hdf5r::H5File$new(h5_path, mode = "r")
  on.exit(f$close_all(), add = TRUE)
  feat <- f[["matrix/features"]]
  idx  <- which(feat[["feature_type"]][] == "CRISPR Guide Capture")
  if (!length(idx)) stop("No CRISPR Guide Capture features found in H5.")
  mat  <- methods::new("dgCMatrix",
                       x   = as.numeric(f[["matrix/data"]][]),
                       i   = f[["matrix/indices"]][],
                       p   = f[["matrix/indptr"]][],
                       Dim = f[["matrix/shape"]][])
  rownames(mat) <- feat[["name"]][]
  colnames(mat) <- f[["matrix/barcodes"]][]
  mat[idx, , drop = FALSE]
}

.auto_threshold <- function(umi_vec, label = "") {
  umi_vec <- umi_vec[umi_vec > 0]
  if (length(umi_vec) < 20) { log_warn("Too few UMIs for auto-threshold (", label, "). Using 5."); return(5) }
  dens   <- stats::density(umi_vec, bw = "SJ", n = 1024)
  dy     <- diff(dens$y)
  valley <- which(dy[-length(dy)] < 0 & dy[-1] >= 0)
  if (!length(valley)) { log_warn("No KDE valley (", label, "). Using median/2."); return(stats::median(umi_vec) / 2) }
  thresh <- dens$x[valley[1] + 1]
  log_info("Auto-threshold [", label, "]: ", round(thresh, 2))
  thresh
}

.get_threshold <- function(CIRI_long, type_filter, fixed_genes, threshold, label) {
  if (threshold != -1) return(threshold)
  umi_vec <- CIRI_long |>
    dplyr::filter(type == type_filter, feature %in% fixed_genes) |>
    dplyr::group_by(cell_barcode) |>
    dplyr::summarise(s = sum(umi), .groups = "drop") |>
    dplyr::pull(s)
  .auto_threshold(umi_vec, label)
}

.fixed_summaries <- function(CIRI_long, a_genes, i_genes, thresh_a, thresh_i) {
  make_s <- function(type_f, fixed_g, label) {
    CIRI_long |> dplyr::filter(type == type_f) |>
      dplyr::mutate(is_fixed = feature %in% fixed_g) |>
      dplyr::group_by(cell_barcode) |>
      dplyr::summarise(total_umi = sum(umi), fixed_umi = sum(umi[is_fixed]),
                       percentage = fixed_umi / total_umi * 100,
                       feature = label, .groups = "drop") |>
      dplyr::filter(total_umi > 0)
  }
  crispra <- make_s("a", a_genes, "CRISPRa")
  crispri <- make_s("i", i_genes, "CRISPRi")
  pos_a   <- dplyr::filter(crispra, fixed_umi >= thresh_a) |> dplyr::pull(cell_barcode)
  pos_i   <- dplyr::filter(crispri, fixed_umi >= thresh_i) |> dplyr::pull(cell_barcode)
  log_info("CRISPRa positive: ", length(pos_a), " | CRISPRi positive: ", length(pos_i))
  list(crispra = crispra, crispri = crispri, pos_a = pos_a, pos_i = pos_i)
}

# ---------------------------------------------------------------------------
# .assign_var_guides: shared helper — assigns one variable guide per cell per type.
# Logic (same as original CIRI pipeline):
#   Among variable guides of a given type, rank by UMI per cell.
#   A cell is assigned if: top guide UMI >= min_umi AND top/second ratio >= min_ratio.
#   Returns a data.frame with cell_barcode and the top guide name for passing cells.
# ---------------------------------------------------------------------------
.assign_var_guides <- function(CIRI_long, fixed_genes, type_filter,
                               min_umi = 10, min_ratio = 5) {
  CIRI_long |>
    dplyr::filter(type == type_filter, !feature %in% fixed_genes) |>
    dplyr::group_by(cell_barcode) |>
    dplyr::arrange(dplyr::desc(umi), .by_group = TRUE) |>
    dplyr::summarise(
      top_guide  = dplyr::first(feature),
      top_umi    = dplyr::first(umi),
      second_umi = ifelse(dplyr::n() >= 2, dplyr::nth(umi, 2), 0),
      ratio      = ifelse(second_umi > 0, top_umi / second_umi, Inf),
      assigned   = top_umi >= min_umi & ratio >= min_ratio,
      .groups    = "drop"
    ) |>
    dplyr::filter(assigned) |>
    dplyr::select(cell_barcode, top_guide)
}

# Strategy 1: assign cells that have EITHER a CRISPRa OR a CRISPRi variable guide
# (or both). Each type is assigned independently.
.assign_single <- function(CIRI_long, a_genes, i_genes, thresh_a, thresh_i) {
  fs    <- .fixed_summaries(CIRI_long, a_genes, i_genes, thresh_a, thresh_i)
  var_a <- .assign_var_guides(CIRI_long, c(a_genes, i_genes), "a") |>
    dplyr::rename(feature_a = top_guide)
  var_i <- .assign_var_guides(CIRI_long, c(a_genes, i_genes), "i") |>
    dplyr::rename(feature_i = top_guide)

  ann <- data.frame(cell_barcode = union(fs$pos_a, fs$pos_i)) |>
    dplyr::left_join(var_a, by = "cell_barcode") |>
    dplyr::left_join(var_i, by = "cell_barcode") |>
    dplyr::mutate(
      feature_a = ifelse(cell_barcode %in% fs$pos_a, tidyr::replace_na(feature_a, "NA"), "NA"),
      feature_i = ifelse(cell_barcode %in% fs$pos_i, tidyr::replace_na(feature_i, "NA"), "NA")
    )
  list(annotation = ann, crispra = fs$crispra, crispri = fs$crispri)
}

# Strategy 2 — CIRI mode: cells must have BOTH a CRISPRa AND a CRISPRi variable guide.
# Each type is assigned independently (top guide >= 10 UMIs, ratio >= 5).
# This is the original "dual" CIRI assignment: one perturbation from each arm.
.assign_dual <- function(CIRI_long, a_genes, i_genes, thresh_a, thresh_i) {
  fs    <- .fixed_summaries(CIRI_long, a_genes, i_genes, thresh_a, thresh_i)
  var_a <- .assign_var_guides(CIRI_long, c(a_genes, i_genes), "a") |>
    dplyr::rename(feature_a = top_guide)
  var_i <- .assign_var_guides(CIRI_long, c(a_genes, i_genes), "i") |>
    dplyr::rename(feature_i = top_guide)

  # CIRI mode: keep only cells that are positive for the fixed guide AND
  # have a valid variable guide assignment for BOTH arms.
  ciri_cells <- intersect(
    intersect(fs$pos_a, var_a$cell_barcode),
    intersect(fs$pos_i, var_i$cell_barcode)
  )
  log_info("CIRI (dual CRISPRa+CRISPRi) cells: ", length(ciri_cells))

  # Also build single-arm annotation for cells positive in only one arm
  only_a <- setdiff(intersect(fs$pos_a, var_a$cell_barcode), fs$pos_i)
  only_i <- setdiff(intersect(fs$pos_i, var_i$cell_barcode), fs$pos_a)
  log_info("CRISPRa-only cells: ", length(only_a),
           " | CRISPRi-only cells: ", length(only_i))

  all_cells <- union(union(ciri_cells, only_a), only_i)

  ann <- data.frame(cell_barcode = all_cells) |>
    dplyr::left_join(var_a, by = "cell_barcode") |>
    dplyr::left_join(var_i, by = "cell_barcode") |>
    dplyr::mutate(
      feature_a = tidyr::replace_na(feature_a, "NA"),
      feature_i = tidyr::replace_na(feature_i, "NA")
    )
  list(annotation = ann, crispra = fs$crispra, crispri = fs$crispri)
}

.save_assignment_plots <- function(CIRI_long, crispra, crispri, plots_dir) {
  plot_df <- CIRI_long |>
    dplyr::group_by(feature, type) |> dplyr::summarise(total_umi = sum(umi), .groups = "drop") |>
    dplyr::arrange(type, feature) |> dplyr::mutate(feature = factor(feature, levels = unique(feature)))
  ggplot2::ggsave(file.path(plots_dir, "total_umixguide.pdf"),
    ggplot2::ggplot(plot_df, ggplot2::aes(x = feature, y = total_umi, fill = type)) +
      ggplot2::geom_bar(stat = "identity") + ggplot2::theme_minimal() +
      ggplot2::labs(title = "Total UMIs per guide") +
      ggplot2::scale_y_continuous(labels = scales::label_number()) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)))

  umis_cell <- CIRI_long |> dplyr::group_by(cell_barcode) |>
    dplyr::summarise(total = sum(umi), .groups = "drop")
  ggplot2::ggsave(file.path(plots_dir, "UMIxcell_hist.pdf"),
    ggplot2::ggplot(umis_cell, ggplot2::aes(x = total)) +
      ggplot2::geom_histogram(binwidth = 1) + ggplot2::labs(title = "UMIs per cell"))

  CRISPRai <- dplyr::bind_rows(crispra, crispri) |>
    dplyr::mutate(variable_umi = total_umi - fixed_umi)
  ggplot2::ggsave(file.path(plots_dir, "fixed_vs_variable_scatter.pdf"),
    ggplot2::ggplot(CRISPRai, ggplot2::aes(x = variable_umi, y = fixed_umi, color = feature)) +
      ggplot2::geom_point(alpha = 0.1, size = 0.5) + ggplot2::theme_minimal() +
      ggplot2::labs(title = "Fixed vs Variable UMI"), width = 12, height = 10)

  sg_dir <- file.path(plots_dir, "single_guide_plots")
  dir.create(sg_dir, showWarnings = FALSE)
  for (guide in unique(CIRI_long$feature)) {
    df <- dplyr::filter(CIRI_long, feature == guide)
    ggplot2::ggsave(file.path(sg_dir, paste0("umi_", guide, ".pdf")),
      ggplot2::ggplot(df, ggplot2::aes(x = umi)) +
        ggplot2::geom_histogram(binwidth = 1) + ggplot2::labs(title = guide))
  }
}
