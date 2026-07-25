# =============================================================================
# CIRI — Guide Name Harmonisation
# =============================================================================

#' Harmonise guide names between guides CSV and CellRanger H5 output
#'
#' CellRanger collapses guide replicate names when building the H5 matrix.
#' For example, \code{ATF7IP_1A} and \code{ATF7IP_1B} in your guides CSV
#' both become \code{ATF7IP_1} in the H5 CRISPR Capture features. If your
#' guides CSV still uses the original A/B names,
#' \code{\link{ciri_step01_assignment}} will fail to match them and return
#' 0 assigned cells.
#'
#' This function reads your original guides CSV and the CellRanger
#' \code{protospacer_calls_per_cell.csv}, figures out the correct name
#' mapping, and writes a harmonised guides CSV ready for Step 01.
#'
#' @section Matching rules (applied in order, first match wins):
#' \enumerate{
#'   \item Strip trailing A or B after \code{_number}:
#'         \code{ATF7IP_1A} → \code{ATF7IP_1}
#'   \item Gene prefix (everything before the first \code{_}):
#'         \code{NTCa_1A} → \code{NTCa}
#'   \item Exact match:
#'         \code{MYOD_1} → \code{MYOD_1}
#' }
#'
#' @param guides_path     Path to the original guides CSV.
#'   Format (no header): \code{feature, type (a/i), fixed (f/v)}.
#' @param protospacer_path Path to \code{protospacer_calls_per_cell.csv}
#'   from CellRanger aggr output.
#' @param out_path        Output path for the harmonised guides CSV.
#'   Default: \code{guides_harmonised.csv} in the same directory as
#'   \code{guides_path}.
#'
#' @return Invisibly returns the harmonised guides \code{data.frame}.
#'   Side effects: writes \code{guides_harmonised.csv} and
#'   \code{name_mapping.csv} to disk.
#'
#' @section Output files:
#' \describe{
#'   \item{\code{guides_harmonised.csv}}{Same 3-column format as the original
#'     guides CSV but with feature names matching the H5 exactly. One row per
#'     unique CellRanger name (A/B replicates deduplicated).}
#'   \item{\code{name_mapping.csv}}{Full mapping table: original name →
#'     CellRanger name, type, fixed, status (matched/unmatched). Use this to
#'     verify the mapping looks correct.}
#' }
#'
#' @examples
#' \dontrun{
#' ciri_harmonise_guides(
#'   guides_path      = "scratch/guides_2.csv",
#'   protospacer_path = "scratch/protospacer_calls_per_cell.csv",
#'   out_path         = "scratch/guides_harmonised.csv"
#' )
#'
#' # Then use the harmonised CSV in Step 01:
#' ciri_step01_assignment(
#'   data_dir = "/path/to/data",
#'   matrix   = "filtered_feature_bc_matrix.h5",
#'   guides   = "scratch/guides_harmonised.csv",
#'   sample   = "AB011",
#'   strategy = 2
#' )
#' }
#'
#' @importFrom utils read.csv write.csv
#' @export
ciri_harmonise_guides <- function(guides_path,
                                   protospacer_path,
                                   out_path = NULL) {

  # ---- Validate inputs -------------------------------------------------------
  if (!file.exists(guides_path))
    stop("guides file not found: ", guides_path, call. = FALSE)
  if (!file.exists(protospacer_path))
    stop("protospacer file not found: ", protospacer_path, call. = FALSE)

  if (is.null(out_path))
    out_path <- file.path(dirname(guides_path), "guides_harmonised.csv")
  mapping_path <- file.path(dirname(out_path), "name_mapping.csv")

  # ---- Banner ----------------------------------------------------------------
  sep <- strrep("=", 72)
  message(sep)
  message("  CIRI — Guide Name Harmonisation")
  message(sprintf("  %s", Sys.time()))
  message(sep)
  message("  Guides CSV   : ", guides_path)
  message("  Protospacer  : ", protospacer_path)
  message("  Output       : ", out_path)
  message("  Name mapping : ", mapping_path)
  message(sep)

  # ---- Load guides -----------------------------------------------------------
  guides <- utils::read.csv(guides_path, header = FALSE,
                             stringsAsFactors = FALSE)
  if (ncol(guides) < 3)
    stop("guides CSV must have 3 columns: feature, type (a/i), fixed (f/v)",
         call. = FALSE)
  names(guides)[1:3] <- c("feature", "type", "fixed")

  # ---- Load protospacer names ------------------------------------------------
  proto <- utils::read.csv(protospacer_path, header = TRUE,
                            stringsAsFactors = FALSE)
  if (!"feature_call" %in% names(proto))
    stop("protospacer CSV must have a 'feature_call' column.\n",
         "Columns found: ", paste(names(proto), collapse = ", "),
         call. = FALSE)

  proto_names <- unique(unlist(strsplit(proto$feature_call, "\\|")))
  proto_names <- proto_names[nchar(proto_names) > 0]

  log_info(length(proto_names), " unique guide names in protospacer CSV")
  log_info(nrow(guides), " guide entries in guides CSV")

  # ---- Name matching ---------------------------------------------------------
  .find_match <- function(name, data_names) {
    # Rule 1: strip trailing A/B after _number  e.g. ATF7IP_1A -> ATF7IP_1
    m <- regmatches(name,
                    regexpr("^.+_\\d+(?=[AB]$)", name, perl = TRUE))
    if (length(m) > 0 && m %in% data_names) return(m)

    # Rule 2: gene prefix (before first _)  e.g. NTCa_1A -> NTCa
    gene <- strsplit(name, "_")[[1]][1]
    if (gene %in% data_names) return(gene)

    # Rule 3: exact match
    if (name %in% data_names) return(name)

    return(NA_character_)
  }

  guides$cellranger_name <- vapply(guides$feature, .find_match,
                                    character(1), data_names = proto_names)

  # ---- Report ----------------------------------------------------------------
  matched   <- guides[!is.na(guides$cellranger_name), ]
  unmatched <- guides[ is.na(guides$cellranger_name), ]

  log_info("Matched   : ", nrow(matched), " / ", nrow(guides), " guides")

  if (nrow(unmatched) > 0) {
    log_warn(nrow(unmatched), " guide(s) not found in protospacer CSV ",
             "(will be excluded):")
    for (nm in unmatched$feature) log_warn("  ", nm)
    message("  Possible reasons:")
    message("    - Guide not detected in this experiment")
    message("    - Name format is unusual — check name_mapping.csv")
  }

  # ---- Build harmonised CSV --------------------------------------------------
  # One row per unique CellRanger name; A/B replicates collapsed.
  harmonised <- list()
  for (cr_name in unique(matched$cellranger_name)) {
    rows <- matched[matched$cellranger_name == cr_name, ]

    if (length(unique(rows$type)) > 1)
      log_warn("Guide '", cr_name, "': conflicting type (",
               paste(unique(rows$type), collapse = "/"), ") — using first")
    if (length(unique(rows$fixed)) > 1)
      log_warn("Guide '", cr_name, "': conflicting fixed (",
               paste(unique(rows$fixed), collapse = "/"), ") — using first")

    harmonised[[cr_name]] <- data.frame(
      feature = cr_name,
      type    = rows$type[1],
      fixed   = rows$fixed[1],
      stringsAsFactors = FALSE
    )
  }

  harmonised_df        <- do.call(rbind, harmonised)
  rownames(harmonised_df) <- NULL

  # ---- Save ------------------------------------------------------------------
  utils::write.csv(harmonised_df, out_path, row.names = FALSE, quote = FALSE)
  log_info("Harmonised guides CSV saved: ", out_path)
  log_info("  ", nrow(harmonised_df), " rows ",
           "(from ", nrow(guides), " original, A/B replicates deduplicated)")

  mapping_df <- data.frame(
    original_name   = guides$feature,
    cellranger_name = guides$cellranger_name,
    type            = guides$type,
    fixed           = guides$fixed,
    status          = ifelse(is.na(guides$cellranger_name),
                             "unmatched", "matched"),
    stringsAsFactors = FALSE
  )
  utils::write.csv(mapping_df, mapping_path, row.names = FALSE, quote = FALSE)
  log_info("Name mapping table saved   : ", mapping_path)

  # ---- Print summary ---------------------------------------------------------
  message("")
  message("=== Harmonised guide annotation ===")
  message(sprintf("  %-18s  %-6s  %-6s", "Feature", "Type", "Fixed"))
  message(sprintf("  %s", strrep("-", 34)))
  for (i in seq_len(nrow(harmonised_df))) {
    message(sprintf("  %-18s  %-6s  %-6s",
                    harmonised_df$feature[i],
                    harmonised_df$type[i],
                    harmonised_df$fixed[i]))
  }
  message("")
  message("Next step — use the harmonised CSV in Step 01:")
  message("  ciri_step01_assignment(")
  message("    ...,")
  message("    guides = \"", out_path, "\",")
  message("    strategy = 2")
  message("  )")

  invisible(harmonised_df)
}
