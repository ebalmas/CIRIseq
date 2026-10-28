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
#' @param strategy    Integer, which cells are retained. \code{1} = every cell
#'                    positive for a fixed guide in either arm is kept, and the
#'                    arm with no valid variable-guide call gets \code{"NA"}.
#'                    \code{2} = a valid variable-guide call is required; dual
#'                    (CRISPRa+CRISPRi) cells are kept together with
#'                    CRISPRa-only and CRISPRi-only cells. Use \code{2} for
#'                    CRISPRa+CRISPRi combinatorial screens.
#'                    This parameter controls the ARMS, not how many variable
#'                    guides are called per arm: for that see
#'                    \code{variable_guides}.
#' @param variable_guides Integer, the variable-guide calling protocol, i.e. how
#'                    many variable guides per arm the library carries. This is
#'                    the \code{guide_strategy} / \code{variable_guides} argument
#'                    of the original \code{CIRI_unified.R}.
#'                    \code{1} = Mode 1, one variable guide per arm: the top
#'                    guide is called when \code{top1 >= min_umi} and
#'                    \code{top1 / top2 >= min_ratio}; \code{feature_a} /
#'                    \code{feature_i} hold a single guide name.
#'                    \code{2} = Mode 2, two variable guides per arm, i.e. the
#'                    A/B pair of the same gene: the top three guides are ranked,
#'                    the pair is called when
#'                    \code{top1 + top2 >= min_umi_pair} and
#'                    \code{(top1 + top2) / top3 >= min_ratio_pair}, and top1 and
#'                    top2 must be the \code{A} and \code{B} guide of the same
#'                    gene and the same number. The two names are returned
#'                    collapsed with a semicolon, e.g. \code{"NTCa_1A;NTCa_1B"}.
#'                    Mode 2 requires the A/B names to be distinct in the input;
#'                    see the note in Details.
#' @param threshold_a CRISPRa fixed-guide UMI threshold. \code{-1} = auto-detect (default).
#' @param threshold_i CRISPRi fixed-guide UMI threshold. \code{-1} = auto-detect (default).
#' @param min_umi     Mode 1 only: minimum UMI count of the top variable guide. Default \code{10}.
#' @param min_ratio   Mode 1 only: minimum \code{top1 / top2} UMI ratio. Default \code{5}.
#' @param min_umi_pair  Mode 2 only: minimum \code{top1 + top2} UMI count. Default \code{4}.
#' @param min_ratio_pair Mode 2 only: minimum \code{(top1 + top2) / top3} UMI ratio. Default \code{10}.
#' @param missing_top3 Mode 2 only, what to do with cells that have exactly two
#'                    variable guides in an arm, so that \code{top3} does not
#'                    exist and the ratio is undefined. \code{"drop"} (default)
#'                    discards them, which is the behaviour of the original
#'                    \code{CIRI_unified.R} and reproduces published results;
#'                    \code{"accept"} treats the ratio as infinite and keeps
#'                    them, subject to the \code{min_umi_pair} threshold. The
#'                    number of cells concerned is reported in the log either way.
#'
#' @param rescue_orphans Logical. With \code{variable_guides = 2}, allow cells of
#'   constructs whose partner guide is undetectable experiment-wide to be called
#'   from the dominant guide alone. Such calls carry ONE guide name; genuine
#'   pairs carry two joined by \code{";"}. They are not validated at this stage -
#'   run \code{ciri_rescue_qc()} after clustering.
#' @param orphan_frac A construct counts as orphan when its weaker guide is
#'   detected in fewer than this fraction of the stronger guide's cells.
#' @param min_cells_construct Constructs whose stronger guide is below this many
#'   cells are treated as dead and never rescued.
#' @param min_umi_rescue,min_ratio_rescue Single-guide criterion applied to
#'   rescued cells.
#' @param rescue_constructs Character vector of guide names to rescue, e.g.
#'   \code{"MEF2C_2B"}. Leave \code{NULL} to rescue every orphan - appropriate
#'   only once each one has passed \code{ciri_rescue_qc()}. On AB014 three of
#'   four orphans failed that check, so an explicit list is the safe default.
#'
#' @return Invisibly returns the path to the output folder.
#'
#' @details
#' \strong{Mode 2 and guide-name harmonisation.} \code{.harmonise_guide_names()}
#' collapses \code{GENE_1A} / \code{GENE_1B} to \code{GENE_1} only when the
#' collapsed name is actually present in the data; when the CellRanger feature
#' reference keeps the two names distinct, the exact names are preserved and
#' Mode 2 can run. Mode 2 therefore needs no change to the harmonisation step,
#' but it does need an input in which the \code{A} and \code{B} guides are
#' separate features: if fewer than two distinct A/B pairs are found in an arm,
#' the function stops with an explanatory error rather than silently returning
#' an empty assignment.
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
                                   strategy        = 1L,
                                   variable_guides = 1L,
                                   threshold_a     = -1,
                                   threshold_i     = -1,
                                   min_umi         = 10,
                                   min_ratio       = 5,
                                   min_umi_pair    = 4,
                                   min_ratio_pair  = 10,
                                   missing_top3    = c("drop", "accept"),
                                   rescue_orphans  = FALSE,
                                   orphan_frac     = 0.10,
                                   min_cells_construct = 50,
                                   min_umi_rescue  = 10,
                                   min_ratio_rescue = 5,
                                   rescue_constructs = NULL) {

  missing_top3 <- match.arg(missing_top3)
  if (rescue_orphans && variable_guides != 2L)
    stop("rescue_orphans applies to variable_guides = 2 only", call. = FALSE)

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
  if (!variable_guides %in% c(1L, 2L))
    stop("variable_guides must be 1 (one variable guide per arm) or 2 (A/B pair per arm)",
         call. = FALSE)

  # ---- load and harmonise guide annotations --------------------------------
  # ciri_harmonise_guides() writes guides_harmonised.csv WITH a header, while the
  # original guides.csv has none: accept both instead of turning the header line
  # into a guide called "feature".
  ann_raw <- read.csv(file.path(data_dir, guides), header = FALSE,
                      stringsAsFactors = FALSE)
  if (nrow(ann_raw) && identical(tolower(trimws(as.character(ann_raw[1, 1]))), "feature")) {
    ann_raw <- ann_raw[-1, , drop = FALSE]
    rownames(ann_raw) <- NULL
    log_info("Header row detected in '", guides, "' and skipped")
  }
  names(ann_raw) <- c("feature", "type", "fixed")

  # Fail before the expensive matrix load if Mode 2 cannot possibly work: the
  # guides CSV must still carry the A/B replicate names. guides_harmonised.csv
  # has them collapsed by design, so Mode 2 needs the pre-harmonisation CSV.
  if (variable_guides == 2L && !any(grepl("[AB]$", ann_raw$feature)))
    stop("variable_guides = 2 (Mode 2) needs a guides CSV whose names still end ",
         "in A/B (e.g. NTCa_1A / NTCa_1B), but '", guides, "' has none. ",
         "guides_harmonised.csv collapses the A/B replicates: pass the original ",
         "guides CSV, and check that the CellRanger feature reference declares ",
         "A and B as separate features.", call. = FALSE)

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
  ann       <- result$ann
  mapping   <- result$mapping

  write.csv(mapping, file.path(out$csv, "guide_name_mapping.csv"), row.names = FALSE)
  log_info("Guide name mapping saved to csv/guide_name_mapping.csv")

  a_genes <- ann[ann$type == "a" & ann$fixed == "f", "feature"]
  i_genes <- ann[ann$type == "i" & ann$fixed == "f", "feature"]
  log_info("Fixed CRISPRa: ", paste(a_genes, collapse = ", "))
  log_info("Fixed CRISPRi: ", paste(i_genes, collapse = ", "))

  write.csv(CIRI_long, file.path(out$csv, "CIRI_long.csv"), row.names = FALSE)

  # ---- thresholds ----------------------------------------------------------
  thresh_a <- .get_threshold(CIRI_long, "a", a_genes, threshold_a, "CRISPRa",
                              plots_dir = out$plots)
  thresh_i <- .get_threshold(CIRI_long, "i", i_genes, threshold_i, "CRISPRi",
                              plots_dir = out$plots)
  log_info("Threshold CRISPRa: ", thresh_a, " | CRISPRi: ", thresh_i)

  # ---- assign --------------------------------------------------------------
  # variable-guide calling protocol: Mode 1 (one guide per arm) or Mode 2 (A/B pair)
  call_var <- if (variable_guides == 1L) {
    log_info("Variable guides: Mode 1, one guide per arm (min_umi = ", min_umi,
             ", min_ratio = ", min_ratio, ")")
    function(cl, fixed_genes, type_filter)
      .assign_var_guides(cl, fixed_genes, type_filter,
                         min_umi = min_umi, min_ratio = min_ratio)
  } else if (!rescue_orphans) {
    log_info("Variable guides: Mode 2, A/B pair per arm (min_umi_pair = ", min_umi_pair,
             ", min_ratio_pair = ", min_ratio_pair, ", missing_top3 = ", missing_top3, ")")
    function(cl, fixed_genes, type_filter)
      .assign_var_guides_pair(cl, fixed_genes, type_filter,
                              min_umi_pair = min_umi_pair,
                              min_ratio_pair = min_ratio_pair,
                              missing_top3 = missing_top3)
  } else {
    log_info("Variable guides: Mode 2 + orphan rescue (orphan_frac = ", orphan_frac,
             ", min_cells = ", min_cells_construct, ", min_umi_rescue = ", min_umi_rescue,
             ", min_ratio_rescue = ", min_ratio_rescue, ")")
    log_warn("Orphan rescue is ON: cells of constructs whose partner guide is ",
             "undetectable experiment-wide are called from a single guide. Such ",
             "calls carry ONE guide name; genuine pairs carry two joined by ';'. ",
             "They are NOT validated here - run ciri_rescue_qc() after clustering.")
    function(cl, fixed_genes, type_filter)
      .assign_var_guides_pair_rescue(cl, fixed_genes, type_filter,
                                     min_umi_pair = min_umi_pair,
                                     min_ratio_pair = min_ratio_pair,
                                     missing_top3 = missing_top3,
                                     min_umi_rescue = min_umi_rescue,
                                     min_ratio_rescue = min_ratio_rescue,
                                     orphan_frac = orphan_frac,
                                     min_cells = min_cells_construct,
                                     rescue_constructs = rescue_constructs)
  }

  assign_result <- if (strategy == 1L)
    .assign_single(CIRI_long, a_genes, i_genes, thresh_a, thresh_i, call_var)
  else
    .assign_dual(CIRI_long, a_genes, i_genes, thresh_a, thresh_i, call_var)

  # ---- plots ---------------------------------------------------------------
  log_info("Saving QC plots ...")
  .save_assignment_plots(CIRI_long, assign_result$crispra, assign_result$crispri, out$plots)

  # ---- stats ---------------------------------------------------------------
  stats_txt <- c(
    paste("Sample            :", sample),
    paste("Strategy (arms)   :", strategy),
    paste("Variable guides   :", variable_guides,
          if (variable_guides == 1L) "(Mode 1: one guide per arm)"
          else if (!rescue_orphans) "(Mode 2: A/B pair per arm)"
          else "(Mode 2 + orphan rescue)"),
    paste("Orphan rescue     :", if (rescue_orphans)
            paste0("ON (orphan_frac = ", orphan_frac, ", min_umi_rescue = ",
                   min_umi_rescue, ", min_ratio_rescue = ", min_ratio_rescue,
                   ", constructs = ",
                   if (is.null(rescue_constructs)) "all orphans" else
                     paste(rescue_constructs, collapse = "/"),
                   "); single-guide calls are those without ';'")
          else "OFF"),
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

.harmonise_guide_names <- function(ann_raw, data_guide_names) {
  find_match <- function(csv_name, data_names) {
    m <- regmatches(csv_name, regexpr("^.+_\\d+(?=[AB]$)", csv_name, perl = TRUE))
    if (length(m) && m %in% data_names) return(m)
    gene <- strsplit(csv_name, "_")[[1]][1]
    if (gene %in% data_names) return(gene)
    if (csv_name %in% data_names) return(csv_name)
    return(NA_character_)
  }

  ann_raw$matched <- vapply(ann_raw$feature, find_match,
                            character(1), data_names = data_guide_names)

  unmatched <- ann_raw$feature[is.na(ann_raw$matched)]
  if (length(unmatched))
    log_warn(length(unmatched), " guide(s) in guides CSV not found in data: ",
             paste(unmatched, collapse = ", "))

  ann_matched <- ann_raw[!is.na(ann_raw$matched), ]
  ann_dedup   <- ann_matched[!duplicated(ann_matched$matched), ]
  ann_out     <- data.frame(feature = ann_dedup$matched,
                            type    = ann_dedup$type,
                            fixed   = ann_dedup$fixed,
                            stringsAsFactors = FALSE)

  mapping <- data.frame(
    guides_csv_name  = ann_raw$feature,
    cellranger_name  = ann_raw$matched,
    type             = ann_raw$type,
    fixed            = ann_raw$fixed,
    stringsAsFactors = FALSE
  )

  list(ann = ann_out, mapping = mapping)
}

.load_from_protospacer <- function(proto_path, ann_raw) {
  proto <- read.csv(proto_path, header = TRUE, stringsAsFactors = FALSE)

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

  data_guide_names <- unique(CIRI_raw$feature)
  harm <- .harmonise_guide_names(ann_raw, data_guide_names)

  CIRI_long <- dplyr::left_join(CIRI_raw, harm$ann, by = "feature") |>
    dplyr::group_by(cell_barcode) |>
    dplyr::mutate(total_umis   = sum(umi),
                  percentage   = umi / total_umis * 100) |>
    dplyr::ungroup()

  list(CIRI_long = CIRI_long, ann = harm$ann, mapping = harm$mapping)
}

.load_from_h5 <- function(h5_path, ann_raw) {
  crispr_mat <- .load_crispr_h5(h5_path)

  CIRI_raw <- as.data.frame(t(as.matrix(crispr_mat))) |>
    tibble::rownames_to_column("cell_barcode") |>
    tidyr::pivot_longer(-cell_barcode, names_to = "feature", values_to = "umi")

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
# Finds the valley in the KDE of per-cell fixed-guide UMI sums.
# Saves a diagnostic PDF to plots_dir if provided.
# ---------------------------------------------------------------------------
.auto_threshold <- function(umi_vec, label = "", plots_dir = NULL) {
  umi_vec <- umi_vec[umi_vec > 0]
  if (length(umi_vec) < 20) {
    log_warn("Too few UMIs for auto-threshold (", label, "). Using 5.")
    return(5)
  }

  dens   <- stats::density(umi_vec, bw = "SJ", n = 1024)
  dy     <- diff(dens$y)
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

  if (!is.null(plots_dir)) {
    df_dens <- data.frame(x = dens$x, y = dens$y)

    p <- ggplot2::ggplot() +
      ggplot2::geom_histogram(
        data    = data.frame(umi = umi_vec),
        mapping = ggplot2::aes(x = umi, y = ggplot2::after_stat(density)),
        bins    = 80, fill = "grey80", colour = "grey60", alpha = 0.6
      ) +
      ggplot2::geom_line(
        data    = df_dens,
        mapping = ggplot2::aes(x = x, y = y),
        colour  = "#2166AC", linewidth = 1
      ) +
      {if (length(all_valleys))
        ggplot2::geom_vline(
          xintercept = dens$x[all_valleys + 1],
          linetype = "dotted", colour = "grey40", linewidth = 0.7
        )
      } +
      ggplot2::geom_vline(
        xintercept = thresh,
        linetype = "dashed", colour = "#D6604D", linewidth = 1.2
      ) +
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
                         paste0("threshold_kde_", gsub(" ", "_", label), ".pdf"))
    ggplot2::ggsave(out_pdf, p, width = 7, height = 5)
    log_info("KDE threshold plot saved: plots/", basename(out_pdf))
  }

  thresh
}

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

# Mode 2: two variable guides per arm, the A/B pair of the same gene.
# Port of the "Dual Guide Logic" branch of CIRI_unified.R (lines 388-470 of
# alessandro-bertero/CIRI_bioinfo), returning the same column contract as
# .assign_var_guides(): cell_barcode + top_guide, where top_guide holds the two
# guide names collapsed with a semicolon, e.g. "NTCa_1A;NTCa_1B".
.assign_var_guides_pair <- function(CIRI_long, fixed_genes, type_filter,
                                    min_umi_pair   = 4,
                                    min_ratio_pair = 10,
                                    missing_top3   = c("drop", "accept")) {
  missing_top3 <- match.arg(missing_top3)
  empty <- data.frame(cell_barcode = character(0), top_guide = character(0),
                      stringsAsFactors = FALSE)

  sec <- CIRI_long |>
    dplyr::filter(type == type_filter, !feature %in% fixed_genes) |>
    dplyr::mutate(gene = stringr::str_extract(feature, "^[^_]+"))
  if (!nrow(sec)) {
    log_warn("Mode 2: no variable guides for arm '", type_filter, "'")
    return(empty)
  }

  # the A/B names must be distinct in the input, otherwise Mode 2 cannot work
  pair_letters <- stringr::str_extract(unique(sec$feature), "[AB]$")
  if (sum(!is.na(pair_letters)) < 2L)
    stop("variable_guides = 2 (Mode 2) requires guide names ending in A and B ",
         "(e.g. NTCa_1A / NTCa_1B), but arm '", type_filter, "' has none. ",
         "The CellRanger feature reference has collapsed the A/B replicates, ",
         "so only variable_guides = 1 can be used on this input.", call. = FALSE)

  # rank the top 3 guides per cell, ties broken in favour of the gene of top1
  top1 <- sec |>
    dplyr::group_by(cell_barcode) |>
    dplyr::arrange(dplyr::desc(umi), feature, .by_group = TRUE) |>
    dplyr::slice(1L) |>
    dplyr::ungroup() |>
    dplyr::select(cell_barcode, top1_gene = gene)

  ranked <- sec |>
    dplyr::left_join(top1, by = "cell_barcode") |>
    dplyr::group_by(cell_barcode) |>
    dplyr::arrange(dplyr::desc(umi), dplyr::desc(gene == top1_gene), feature,
                   .by_group = TRUE) |>
    dplyr::mutate(rank = dplyr::row_number()) |>
    dplyr::ungroup() |>
    dplyr::filter(rank <= 3L) |>
    dplyr::select(cell_barcode, rank, umi, feature) |>
    tidyr::pivot_wider(names_from = rank, values_from = c(umi, feature),
                       names_prefix = "top")

  for (nm in c("umi_top1", "umi_top2", "umi_top3", "feature_top1", "feature_top2"))
    if (!nm %in% names(ranked))
      ranked[[nm]] <- if (startsWith(nm, "umi")) NA_integer_ else NA_character_

  n_no_top3 <- sum(!is.na(ranked$umi_top2) & is.na(ranked$umi_top3))
  log_info("Mode 2, arm '", type_filter, "': ", n_no_top3,
           " cell(s) with exactly two variable guides (no top3) -> ", missing_top3)

  ranked <- ranked |>
    dplyr::mutate(
      ratio = dplyr::case_when(
        is.na(umi_top2) | umi_top2 <= 0                     ~ Inf,
        is.na(umi_top3) & missing_top3 == "accept"          ~ Inf,
        TRUE                                                ~ (umi_top1 + umi_top2) / umi_top3
      ),
      pass_filter = !is.na(ratio) & (umi_top1 + umi_top2) >= min_umi_pair &
                    ratio >= min_ratio_pair
    ) |>
    dplyr::filter(pass_filter)

  if (!nrow(ranked)) {
    log_warn("Mode 2: no cell passed the pair filter for arm '", type_filter, "'")
    return(empty)
  }

  # top1 and top2 must be the A and B guide of the same gene and the same number
  paired <- ranked |>
    dplyr::mutate(
      gene1   = stringr::str_extract(feature_top1, "^[^_]+"),
      gene2   = stringr::str_extract(feature_top2, "^[^_]+"),
      num1    = stringr::str_extract(feature_top1, "(?<=_)\\d+"),
      num2    = stringr::str_extract(feature_top2, "(?<=_)\\d+"),
      letter1 = stringr::str_extract(feature_top1, "[A-Za-z]$"),
      letter2 = stringr::str_extract(feature_top2, "[A-Za-z]$")
    ) |>
    dplyr::filter(
      (is.na(feature_top1) & !is.na(feature_top2)) |
      (!is.na(feature_top1) & is.na(feature_top2)) |
      (gene1 == gene2 & num1 == num2 &
         ((letter1 == "A" & letter2 == "B") | (letter1 == "B" & letter2 == "A")))
    )

  log_info("Mode 2, arm '", type_filter, "': ", nrow(paired), " of ", nrow(ranked),
           " cell(s) carry a matching A/B pair")

  paired |>
    dplyr::mutate(top_guide = vapply(
      seq_len(dplyr::n()),
      function(k) {
        f <- sort(stats::na.omit(c(feature_top1[k], feature_top2[k])))
        if (!length(f)) NA_character_ else paste(f, collapse = ";")
      }, character(1))) |>
    dplyr::filter(!is.na(top_guide)) |>
    dplyr::select(cell_barcode, top_guide)
}

.assign_single <- function(CIRI_long, a_genes, i_genes, thresh_a, thresh_i,
                           call_var = .assign_var_guides) {
  fs    <- .fixed_summaries(CIRI_long, a_genes, i_genes, thresh_a, thresh_i)
  var_a <- call_var(CIRI_long, c(a_genes, i_genes), "a") |>
    dplyr::rename(feature_a = top_guide)
  var_i <- call_var(CIRI_long, c(a_genes, i_genes), "i") |>
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

.assign_dual <- function(CIRI_long, a_genes, i_genes, thresh_a, thresh_i,
                         call_var = .assign_var_guides) {
  fs    <- .fixed_summaries(CIRI_long, a_genes, i_genes, thresh_a, thresh_i)
  var_a <- call_var(CIRI_long, c(a_genes, i_genes), "a") |>
    dplyr::rename(feature_a = top_guide)
  var_i <- call_var(CIRI_long, c(a_genes, i_genes), "i") |>
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
