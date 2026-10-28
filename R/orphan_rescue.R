# =============================================================================
# CIRI - Rescue of dual-guide constructs with an undetectable partner guide
# =============================================================================
#
# Added in 0.4.0. See NEWS.md and the validation on AB014.

## Rescue of dual-guide constructs whose partner guide is systematically
## undetectable (silencing, poor capture, or a barcode that was never read out).
##
## Rationale. The Mode 2 pairing filter exists to confirm vector identity: two
## guides of one construct must be seen together, which ambient guide RNA cannot
## easily fake. When one guide of a construct is unobservable ACROSS THE WHOLE
## EXPERIMENT, that confirmation is unavailable for every cell carrying it, and
## requiring it discards real cells (in AB014: all MEF2C cells, and most NTCa).
##
## The distinction this code rests on:
##   - partner missing IN A CELL, for a construct whose partner is normally seen
##     -> informative, the cell probably carries ambient RNA. NOT rescued.
##   - partner missing FOR THE WHOLE CONSTRUCT
##     -> uninformative, the pairing test cannot be applied at all. Rescued,
##        using the single-guide (Mode 1) criterion, and flagged as such.
##
## Rescued calls carry ONE guide name; strict pair calls carry two joined by
## ";". No schema change, and `grepl(";", feature_a)` separates them downstream.
## Gene-level columns (gene_a / gene_i / gene_comb) are unaffected, because the
## gene is parsed from the name either way.

## ---------------------------------------------------------------- step 1 ----
## Classify every A/B construct of one arm from the data.
##   orphan_frac : the weaker guide must be detected in fewer than this fraction
##                 of the stronger guide's cells for the construct to count as
##                 orphan (default 0.10)
##   min_cells   : a construct whose stronger guide is below this is "dead" -
##                 nothing to rescue, and rescuing noise would be worse
.classify_constructs <- function(CIRI_long, fixed_genes, type_filter,
                                 orphan_frac = 0.10, min_cells = 50) {
  sec <- CIRI_long |>
    dplyr::filter(type == type_filter,
                  !feature %in% fixed_genes,
                  grepl("_[0-9]+[AB]$", feature),
                  umi > 0)
  if (!nrow(sec)) return(tibble::tibble())

  det <- sec |>
    dplyr::mutate(construct = sub("[AB]$", "", feature),
                  letter    = sub("^.*_[0-9]+", "", feature)) |>
    dplyr::count(construct, letter, name = "cells") |>
    tidyr::complete(construct, letter = c("A", "B"), fill = list(cells = 0L))

  cells_in <- sec |>
    dplyr::mutate(construct = sub("[AB]$", "", feature)) |>
    dplyr::distinct(construct, cell_barcode, feature) |>
    dplyr::count(construct, cell_barcode, name = "n_guides") |>
    dplyr::filter(n_guides == 2) |>
    dplyr::count(construct, name = "cells_both")

  det |>
    tidyr::pivot_wider(names_from = letter, values_from = cells,
                       values_fill = 0, names_prefix = "cells_") |>
    dplyr::left_join(cells_in, by = "construct") |>
    dplyr::mutate(
      cells_both = tidyr::replace_na(cells_both, 0L),
      strong     = pmax(cells_A, cells_B),
      weak       = pmin(cells_A, cells_B),
      weak_frac  = ifelse(strong > 0, weak / strong, NA_real_),
      dominant   = ifelse(cells_A >= cells_B, paste0(construct, "A"),
                                              paste0(construct, "B")),
      class = dplyr::case_when(
        strong < min_cells                ~ "dead",
        weak_frac < orphan_frac           ~ "orphan",
        TRUE                              ~ "paired")) |>
    dplyr::arrange(class, dplyr::desc(strong))
}

## ---------------------------------------------------------------- step 2 ----
## Two-tier variable-guide caller: strict pairs first, then single-guide rescue
## restricted to the dominant guide of orphan constructs.
## Drop-in replacement for the `call_var` callback of ciri_step01_assignment():
##   same (CIRI_long, fixed_genes, type_filter) signature,
##   same cell_barcode + top_guide contract.
.assign_var_guides_pair_rescue <- function(CIRI_long, fixed_genes, type_filter,
                                           min_umi_pair    = 4,
                                           min_ratio_pair  = 10,
                                           missing_top3    = "drop",
                                           min_umi_rescue  = 10,
                                           min_ratio_rescue = 5,
                                           orphan_frac     = 0.10,
                                           min_cells       = 50,
                                           rescue_constructs = NULL,
                                           report          = TRUE) {

  ## tier 1 - the published Mode 2 logic, untouched
  pairs <- .assign_var_guides_pair(CIRI_long, fixed_genes, type_filter,
                                   min_umi_pair   = min_umi_pair,
                                   min_ratio_pair = min_ratio_pair,
                                   missing_top3   = missing_top3)

  cls <- .classify_constructs(CIRI_long, fixed_genes, type_filter,
                              orphan_frac = orphan_frac, min_cells = min_cells)
  orphans <- cls$dominant[cls$class == "orphan"]

  ## Being an orphan makes a construct ELIGIBLE for rescue; it does not make the
  ## rescue correct. A guide that is merely very abundant wins the single-guide
  ## ratio test in cells that never carried it, and the recovered cells then
  ## track the background composition instead of the construct's real phenotype.
  ## Pass rescue_constructs to restrict the rescue to the constructs that passed
  ## that check (see ciri_validate_rescue()); leaving it NULL rescues every orphan and
  ## is only appropriate when the check has been done for all of them.
  if (!is.null(rescue_constructs)) {
    asked   <- setdiff(rescue_constructs, orphans)
    if (length(asked))
      warning("rescue_constructs not classified as orphan in arm '", type_filter,
              "', ignored: ", paste(asked, collapse = ", "), call. = FALSE)
    dropped <- setdiff(orphans, rescue_constructs)
    if (length(dropped) && report)
      message("[rescue] arm '", type_filter, "': orphan but NOT rescued by request: ",
              paste(sort(dropped), collapse = ", "))
    orphans <- intersect(orphans, rescue_constructs)
  }

  if (!length(orphans)) {
    if (report) message("[rescue] arm '", type_filter,
                        "': no orphan construct, nothing to rescue")
    return(pairs)
  }

  ## tier 2 - single-guide criterion, only for cells tier 1 could not call, and
  ## only when the dominant guide of an orphan construct is the top guide
  sec <- CIRI_long |>
    dplyr::filter(type == type_filter, !feature %in% fixed_genes,
                  !cell_barcode %in% pairs$cell_barcode)

  resc <- sec |>
    dplyr::group_by(cell_barcode) |>
    dplyr::arrange(dplyr::desc(umi), feature, .by_group = TRUE) |>
    dplyr::summarise(top1     = dplyr::first(feature),
                     umi_top1 = dplyr::first(umi),
                     umi_top2 = dplyr::nth(umi, 2, default = 0),
                     .groups  = "drop") |>
    dplyr::filter(top1 %in% orphans,
                  umi_top1 >= min_umi_rescue,
                  umi_top2 == 0 | umi_top1 / umi_top2 >= min_ratio_rescue) |>
    dplyr::transmute(cell_barcode, top_guide = top1)

  if (report) {
    message("[rescue] arm '", type_filter, "': orphan construct(s) ",
            paste(sort(orphans), collapse = ", "))
    message("[rescue] arm '", type_filter, "': ", nrow(pairs),
            " pair call(s) + ", nrow(resc), " single-guide rescue(s)")
  }

  dplyr::bind_rows(pairs, resc)
}

## --------------------------------------------------------------- step 3 ----
## Does a rescued population behave like the construct, or like the background?
##
## This is the check that decides whether a rescue is legitimate. For each
## construct it compares, WITHIN each sample, the phenotype of the rescued cells
## against (a) the pair-confirmed cells of the same construct, where they exist,
## and (b) the sample's overall composition. A rescued set that tracks the
## background is ambient contamination by an abundant barcode, not a recovered
## population, and the construct must be left out of rescue_constructs.
##
##   annotation : annotation_data.csv read with na.strings = character(0)
##   cell_meta  : data.frame(cell, sample, phenotype) where phenotype is logical
##                (e.g. cluster %in% c("3","8"))
#' Does a rescued population behave like the construct, or like the background?
#'
#' The check that decides whether a rescue is legitimate. Being an orphan makes a
#' construct eligible for rescue; only this comparison makes it correct.
#'
#' @param annotation \code{annotation_data.csv} read with
#'   \code{na.strings = character(0)} so the literal \code{"NA"} survives.
#' @param cell_meta \code{data.frame(cell, sample, phenotype)} with
#'   \code{phenotype} logical, e.g. \code{cluster \%in\% c("3","8")}.
#' @param arms Arms to validate.
#' @return A \code{data.frame} with pair, rescued and background percentages
#'   per construct and sample, plus a verdict.
#' @export
ciri_validate_rescue <- function(annotation, cell_meta, arms = c("a", "i")) {
  calls <- dplyr::bind_rows(lapply(arms, function(ar) {
    col <- paste0("feature_", ar)
    annotation[annotation[[col]] != "NA", , drop = FALSE] |>
      dplyr::transmute(cell = cell_barcode, arm = ar, call = .data[[col]])
  })) |>
    dplyr::mutate(kind = dplyr::if_else(grepl(";", call), "pair", "rescued"),
                  # key on the construct, not the guide: a pair is "GENE_nA;GENE_nB"
                  # and a rescue is whichever of the two is dominant, so both must
                  # collapse to "GENE_n" or they never line up.
                  construct = sub("[AB]$", "", sub(";.*$", "", call))) |>
    dplyr::inner_join(cell_meta, by = "cell")

  bg <- cell_meta |>
    dplyr::group_by(sample) |>
    dplyr::summarise(n_background = dplyr::n(),
                     pct_background = 100 * mean(phenotype), .groups = "drop")

  calls |>
    dplyr::group_by(sample, arm, construct, kind) |>
    dplyr::summarise(n = dplyr::n(), pct = 100 * mean(phenotype), .groups = "drop") |>
    dplyr::mutate(kind = factor(kind, levels = c("pair", "rescued"))) |>
    tidyr::pivot_wider(names_from = kind, values_from = c(n, pct),
                       values_fill = NA, names_expand = TRUE) |>
    dplyr::left_join(bg, by = "sample") |>
    dplyr::mutate(dplyr::across(dplyr::starts_with("pct"), ~ round(.x, 1)),
                  # The only positive evidence is resembling the construct's own
                  # pair-confirmed cells. Merely differing from the background is
                  # not evidence of anything - NTCa_1A in AB014 differs from the
                  # background and is still contamination. Where no pairs exist
                  # the verdict is deliberately left to the reader.
                  verdict = dplyr::case_when(
                    is.na(pct_rescued)                      ~ "no rescued cells",
                    is.na(pct_pair)                         ~ "no pair reference - judge by enrichment",
                    abs(pct_rescued - pct_pair) <= 10       ~ "consistent with its own pairs",
                    TRUE ~ "inconsistent with its own pairs - DO NOT rescue"),
                  fold_vs_background = round(pct_rescued / pmax(pct_background, 0.1), 2))
}

## --------------------------------------------------------------- step 4 ----
## The validation as a figure. Note this CANNOT live inside step 01: it needs a
## per-cell phenotype (cluster, score, pseudotime) that only exists after
## clustering. Run it after step 02+, and before trusting any rescued call.
#' Plot the rescue validation
#'
#' @param validation Output of \code{ciri_validate_rescue()}.
#' @param phenotype_label Axis label describing the phenotype.
#' @param file Optional path; when given the figure is written there.
#' @param width,height Figure size in inches.
#' @return A \code{ggplot} object.
#' @export
ciri_plot_rescue_validation <- function(validation, phenotype_label = "cells in mature myotubes (%)",
                                   file = NULL, width = 9.4, height = 3.8) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  d <- validation |>
    dplyr::filter(!is.na(pct_rescued)) |>
    tidyr::pivot_longer(c(pct_pair, pct_rescued), names_to = "kind", values_to = "pct") |>
    dplyr::mutate(
      n = ifelse(kind == "pct_pair", n_pair, n_rescued),
      kind = factor(ifelse(kind == "pct_pair", "pair-confirmed", "single-guide rescued"),
                    levels = c("pair-confirmed", "single-guide rescued")),
      label = ifelse(is.na(pct), "no pairs exist", paste0("n=", n)),
      pct_plot = ifelse(is.na(pct), 0, pct))

  p <- ggplot2::ggplot(d, ggplot2::aes(x = pct_plot, y = construct, fill = kind)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.75),
                      width = 0.7, colour = NA) +
    ggplot2::geom_vline(data = dplyr::distinct(d, sample, pct_background),
                        ggplot2::aes(xintercept = pct_background),
                        linetype = "dashed", colour = "grey55", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = label, x = 1.5),
                       position = ggplot2::position_dodge(width = 0.75),
                       hjust = 0, size = 2.2, colour = "grey20") +
    ggplot2::facet_wrap(~ sample, nrow = 1) +
    ggplot2::scale_fill_manual(values = c("pair-confirmed" = "#1F6FB4",
                                          "single-guide rescued" = "#D1495B"),
                               name = NULL) +
    ggplot2::labs(x = phenotype_label, y = NULL,
                  subtitle = paste("Dashed line = sample background.",
                                   "A rescued population is credible only when it",
                                   "matches the construct's own pairs.")) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(legend.position = "top", panel.grid.major.y = ggplot2::element_blank(),
                   strip.background = ggplot2::element_rect(fill = "grey95", colour = "grey70"),
                   strip.text = ggplot2::element_text(face = "bold"))

  if (!is.null(file)) ggplot2::ggsave(file, p, width = width, height = height)
  p
}

## --------------------------------------------------------------- step 5 -----
## Human-readable audit of what the rescue did, to be written next to
## annotation_data.csv so the relaxation is never invisible.
#' Classify A/B constructs as paired, orphan or dead
#'
#' A construct is \emph{orphan} when its weaker guide is detected in fewer than
#' \code{orphan_frac} of the stronger guide's cells, i.e. the pairing test of
#' Mode 2 cannot be applied to any cell carrying it.
#'
#' @param CIRI_long Long guide table as built by \code{ciri_step01_assignment()}.
#' @param fixed_genes Character vector of fixed-guide feature names.
#' @param arms Arms to report on.
#' @param ... Passed to the classifier (\code{orphan_frac}, \code{min_cells}).
#' @return A \code{data.frame}, one row per construct per arm.
#' @export
ciri_rescue_report <- function(CIRI_long, fixed_genes, arms = c("a", "i"), ...) {
  dplyr::bind_rows(lapply(arms, function(arm) {
    .classify_constructs(CIRI_long, fixed_genes, arm, ...) |>
      dplyr::mutate(arm = arm, .before = 1)
  }))
}


## --------------------------------------------------------------- step 6 ----
#' Rescue QC: validate, plot, and fail if a rescued construct does not hold up
#'
#' Run this after clustering, not inside \code{ciri_step01_assignment()}: the
#' validation needs a per-cell phenotype that does not exist at assignment time.
#'
#' @param annotation_path Path to \code{annotation_data.csv} from step 01.
#' @param cell_meta \code{data.frame(cell, sample, phenotype)}; \code{phenotype}
#'   logical. Any per-cell quantity that separates a real perturbation
#'   population from the background will do - it is not restricted to clusters.
#' @param out_dir Directory for \code{rescue_validation.csv} and
#'   \code{rescue_validation.pdf}.
#' @param rescue_constructs The constructs that were rescued, i.e. what was
#'   passed to \code{ciri_step01_assignment()}. Each is checked, and the
#'   function stops if one is inconsistent with its own pair-confirmed cells.
#' @param phenotype_label Axis label for the figure.
#' @param strict If \code{TRUE} (default) a failing construct is an error; set
#'   \code{FALSE} to downgrade to a warning.
#' @return Invisibly, the validation \code{data.frame}.
#' @export
ciri_rescue_qc <- function(annotation_path, cell_meta, out_dir,
                           rescue_constructs = NULL,
                           phenotype_label = "cells with the phenotype (%)",
                           strict = TRUE) {
  stopifnot(file.exists(annotation_path))
  if (!all(c("cell", "sample", "phenotype") %in% names(cell_meta)))
    stop("cell_meta needs columns cell, sample, phenotype", call. = FALSE)
  if (!is.logical(cell_meta$phenotype))
    stop("cell_meta$phenotype must be logical", call. = FALSE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  annotation <- utils::read.csv(annotation_path, stringsAsFactors = FALSE,
                                na.strings = character(0))
  val <- ciri_validate_rescue(annotation, cell_meta)
  utils::write.csv(val, file.path(out_dir, "rescue_validation.csv"), row.names = FALSE)

  if (any(!is.na(val$pct_rescued))) {
    ciri_plot_rescue_validation(val, phenotype_label,
                                file = file.path(out_dir, "rescue_validation.pdf"))
  } else {
    log_info("No single-guide rescued calls in this annotation, nothing to plot")
  }

  if (!is.null(rescue_constructs)) {
    want <- sub("[AB]$", "", rescue_constructs)
    bad  <- val[val$construct %in% want &
                  grepl("DO NOT rescue", val$verdict), , drop = FALSE]
    if (nrow(bad)) {
      txt <- paste0(bad$construct, " in ", bad$sample, " (rescued ",
                    bad$pct_rescued, "% vs pairs ", bad$pct_pair, "%)",
                    collapse = "; ")
      m <- paste0("rescued construct(s) inconsistent with their own pairs: ", txt)
      if (strict) stop(m, call. = FALSE) else log_warn(m)
    } else {
      log_info("No rescued constructs were found inconsistent with their own pairs; ",
               "constructs without a pair reference still require manual assessment")
    }
  }
  invisible(val)
}
