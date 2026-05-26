# =============================================================================
# CIRI — Step 01: Perturbation Assignment
# =============================================================================

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
#' @param strategy    Integer. \code{1} = single variable guide (default),
#'                    \code{2} = dual variable guides.
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
  thresh_a <- .get_threshold(CIRI_long, "a", a_genes, threshold_a, "CRISPRa")
  thresh_i <- .get_threshold(CIRI_long, "i", i_genes, threshold_i, "CRISPRi")
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
            file.path(out$to_scratch, "annotation_data.csv"), row.names = FALSE)

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

.assign_single <- function(CIRI_long, a_genes, i_genes, thresh_a, thresh_i) {
  fs  <- .fixed_summaries(CIRI_long, a_genes, i_genes, thresh_a, thresh_i)
  var <- CIRI_long |> dplyr::filter(!feature %in% c(a_genes, i_genes)) |>
    dplyr::group_by(cell_barcode, type) |> dplyr::arrange(dplyr::desc(umi), .by_group = TRUE) |>
    dplyr::summarise(top_guide = dplyr::first(feature), top_umi = dplyr::first(umi),
                     second_umi = ifelse(dplyr::n() >= 2, dplyr::nth(umi, 2), 0),
                     ratio = ifelse(second_umi > 0, top_umi / second_umi, Inf),
                     assigned = top_umi >= 10 & ratio >= 5, .groups = "drop")
  var_a <- dplyr::filter(var, type == "a", assigned) |> dplyr::select(cell_barcode, feature_a = top_guide)
  var_i <- dplyr::filter(var, type == "i", assigned) |> dplyr::select(cell_barcode, feature_i = top_guide)
  ann <- data.frame(cell_barcode = union(fs$pos_a, fs$pos_i)) |>
    dplyr::left_join(var_a, by = "cell_barcode") |>
    dplyr::left_join(var_i, by = "cell_barcode") |>
    dplyr::mutate(
      feature_a = ifelse(cell_barcode %in% fs$pos_a, tidyr::replace_na(feature_a, "NA"), "NA"),
      feature_i = ifelse(cell_barcode %in% fs$pos_i, tidyr::replace_na(feature_i, "NA"), "NA"))
  list(annotation = ann, crispra = fs$crispra, crispri = fs$crispri)
}

.assign_dual <- function(CIRI_long, a_genes, i_genes, thresh_a, thresh_i) {
  fs  <- .fixed_summaries(CIRI_long, a_genes, i_genes, thresh_a, thresh_i)
  var <- CIRI_long |> dplyr::filter(!feature %in% c(a_genes, i_genes)) |>
    dplyr::group_by(cell_barcode, type) |> dplyr::arrange(dplyr::desc(umi), .by_group = TRUE) |>
    dplyr::summarise(
      top_guide    = dplyr::first(feature),
      second_guide = ifelse(dplyr::n() >= 2, dplyr::nth(feature, 2), NA_character_),
      top_umi      = dplyr::first(umi),
      second_umi   = ifelse(dplyr::n() >= 2, dplyr::nth(umi, 2), 0),
      third_umi    = ifelse(dplyr::n() >= 3, dplyr::nth(umi, 3), 0),
      top2_sum     = top_umi + second_umi,
      ratio        = ifelse(third_umi > 0, top2_sum / third_umi, Inf),
      same_gene    = !is.na(second_guide) &
        stringr::str_extract(top_guide, "^[^_]+") == stringr::str_extract(second_guide, "^[^_]+"),
      assigned     = top2_sum >= 4 & ratio >= 10 & same_gene,
      combined     = ifelse(assigned, paste(top_guide, second_guide, sep = ";"), NA_character_),
      .groups = "drop")
  var_a <- dplyr::filter(var, type == "a", assigned) |> dplyr::select(cell_barcode, feature_a = combined)
  var_i <- dplyr::filter(var, type == "i", assigned) |> dplyr::select(cell_barcode, feature_i = combined)
  ann <- data.frame(cell_barcode = union(fs$pos_a, fs$pos_i)) |>
    dplyr::left_join(var_a, by = "cell_barcode") |>
    dplyr::left_join(var_i, by = "cell_barcode") |>
    dplyr::mutate(feature_a = tidyr::replace_na(feature_a, "NA"),
                  feature_i = tidyr::replace_na(feature_i, "NA"))
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
