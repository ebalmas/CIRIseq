# =============================================================================
# CIRI — Guide Name Harmonisation
# =============================================================================

#' Harmonise guide names between guides CSV and CellRanger H5 output
#'
#' Some CellRanger references collapse guide replicate names in the H5 matrix.
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
#'   \item Exact match against the matrix feature names.
#'   \item Feature-reference lookup when matching against feature ids.
#'   \item Strip trailing A or B after \code{_number}, provided the result
#'         does not identify a different guide.
#'   \item Gene prefix (everything before the first \code{_}).
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
#' @param feature_reference_path Optional path to the CellRanger
#'   \code{feature_reference.csv}. Supplying it lets the matching use the
#'   \code{id} / \code{name} correspondence instead of guessing from the
#'   protospacer names, and enables the guard that refuses to map a guide onto
#'   the \code{id} of a different guide. Strongly recommended whenever
#'   \code{id} and \code{name} are not the same column.
#' @param matrix_path Optional path to \code{filtered_feature_bc_matrix.h5}.
#'   When given, guide names are matched against \code{features/name} from the
#'   matrix itself, which is what the rest of the pipeline consumes. This is the
#'   most reliable option and takes precedence over
#'   \code{feature_reference_path}. Requires \pkg{hdf5r}.
#' @param allow_collapse Controls what happens when A/B replicate names would be
#'   collapsed. \code{NA} (default) decides from the data: an error if the matrix
#'   resolves A and B as separate features, a warning otherwise. \code{TRUE}
#'   always collapses, \code{FALSE} always errors.
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
                                   out_path = NULL,
                                   feature_reference_path = NULL,
                                   matrix_path = NULL,
                                   allow_collapse = NA) {

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

  # ---- Which names does the MATRIX actually use? -----------------------------
  # CellRanger writes two different columns of the feature reference into two
  # different outputs: features/name into the H5, and id into
  # crispr_analysis/protospacer_calls_per_cell.csv. When id and name differ -
  # e.g. id "ATF7IP_1" for name "ATF7IP_2A", which is how a 4-guides-per-gene
  # reference numbers its rows - matching guide names against protospacer ids
  # silently maps a guide onto a DIFFERENT guide's feature. The rest of the
  # pipeline consumes the H5, so the H5's names are authoritative.
  id2name <- NULL
  if (!is.null(feature_reference_path)) {
    if (!file.exists(feature_reference_path))
      stop("feature reference not found: ", feature_reference_path, call. = FALSE)
    fref <- utils::read.csv(feature_reference_path, stringsAsFactors = FALSE)
    if (!all(c("id", "name") %in% names(fref)))
      stop("feature reference must have 'id' and 'name' columns", call. = FALSE)
    id2name <- stats::setNames(fref$name, fref$id)
    log_info("Feature reference: ", nrow(fref), " rows, id differs from name in ",
             sum(fref$id != fref$name))
  }

  data_names <- proto_names
  data_source <- "protospacer CSV"

  if (!is.null(matrix_path)) {
    if (!requireNamespace("hdf5r", quietly = TRUE))
      stop("matrix_path needs the hdf5r package", call. = FALSE)
    h5 <- hdf5r::H5File$new(matrix_path, mode = "r")
    on.exit(h5$close_all(), add = TRUE)
    ftype <- h5[["matrix/features/feature_type"]][]
    fname <- h5[["matrix/features/name"]][]
    data_names  <- unique(fname[ftype == "CRISPR Guide Capture"])
    data_source <- "H5 features/name"
  } else if (!is.null(id2name) && any(proto_names %in% names(id2name))) {
    translated <- unname(id2name[proto_names])
    n_tr <- sum(!is.na(translated) & translated != proto_names)
    data_names  <- ifelse(is.na(translated), proto_names, translated)
    data_source <- paste0("protospacer CSV, ", n_tr, " id(s) translated to name")
  }
  log_info("Matching guides against: ", data_source, " (", length(data_names), " names)")

  # ---- Name matching ---------------------------------------------------------

  guides$cellranger_name <- vapply(guides$feature, .find_match,
                                    character(1), data_names = data_names,
                                    id2name = id2name)

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

  # ---- Would this destroy a dual-guide library? ------------------------------
  # Collapsing A/B is correct only when the matrix has ONE feature per construct.
  # On a matrix that resolves A and B separately it removes exactly the
  # information Mode 2 (variable_guides = 2) needs, and the failure downstream is
  # an unhelpful "needs a guides CSV whose names still end in A/B".
  ab_in_guides <- sum(grepl("_[0-9]+[AB]$", guides$feature))
  ab_kept      <- sum(grepl("_[0-9]+[AB]$", matched$cellranger_name))
  if (ab_in_guides > 0 && ab_kept < ab_in_guides) {
    msg <- paste0(ab_in_guides - ab_kept, " of ", ab_in_guides,
                  " A/B guide names would be collapsed")
    if (isTRUE(allow_collapse)) {
      log_warn(msg, " - proceeding because allow_collapse = TRUE")
    } else if (identical(allow_collapse, FALSE)) {
      stop(msg, ", and allow_collapse = FALSE.", call. = FALSE)
    } else if (any(grepl("_[0-9]+[AB]$", data_names))) {
      stop(msg, ", but the matrix resolves A and B as separate features ",
           "(e.g. ", paste(utils::head(grep("_[0-9]+[AB]$", data_names, value = TRUE), 2),
                           collapse = ", "),
           "). Harmonisation is not needed for this input and would make ",
           "variable_guides = 2 impossible: pass the guides CSV to ",
           "ciri_step01_assignment() unchanged. Set allow_collapse = TRUE to ",
           "override.", call. = FALSE)
    } else {
      log_warn(msg, " - the matrix has no A/B names, so this is the collapsed ",
               "case harmonisation exists for. variable_guides = 2 will not be ",
               "available on this input.")
    }
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

.find_match <- function(name, data_names, id2name = NULL) {
  # Rule 1: exact match wins. A guide name that is present in the matrix is
  # never a candidate for collapsing - this used to be tried last.
  if (name %in% data_names) return(name)

  # Rule 2: reverse lookup through the feature reference. When the data names
  # are ids, the id carrying THIS guide's name is the right answer, and no
  # string surgery can find it: id "ATF7IP_1" is name "ATF7IP_2A".
  if (!is.null(id2name)) {
    hit <- names(id2name)[!is.na(id2name) & id2name == name]
    hit <- hit[hit %in% data_names]
    if (length(hit) == 1) return(hit)
  }

  # Rule 3: strip trailing A/B after _number  e.g. ATF7IP_1A -> ATF7IP_1, for
  # a reference that really did collapse the two replicates into one feature.
  # Refused when the candidate is the id of a DIFFERENT A/B guide - which is
  # what makes the plain string rule unsafe. An id whose name is the collapsed
  # name itself is the legitimate case and still passes.
  m <- regmatches(name, regexpr("^.+_\\d+(?=[AB]$)", name, perl = TRUE))
  if (length(m) > 0 && m %in% data_names) {
    clash <- !is.null(id2name) && m %in% names(id2name) &&
      !is.na(id2name[[m]]) && id2name[[m]] != name &&
      grepl("_[0-9]+[AB]$", id2name[[m]])
    if (clash) {
      log_warn("Refusing to map '", name, "' to '", m, "': '", m,
               "' is the id of a different guide ('", id2name[[m]], "')")
    } else {
      return(m)
    }
  }

  # Rule 4: gene prefix (before first _)  e.g. NTCa_1A -> NTCa
  gene <- strsplit(name, "_")[[1]][1]
  if (gene %in% data_names) return(gene)

  return(NA_character_)
}
