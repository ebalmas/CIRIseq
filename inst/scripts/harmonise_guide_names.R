#!/usr/bin/env Rscript
# =============================================================================
# harmonise_guide_names.R
# =============================================================================
# PURPOSE
# -------
# CellRanger collapses guide replicate names when building the H5 matrix.
# For example, "ATF7IP_1A" and "ATF7IP_1B" in your guides CSV both become
# "ATF7IP_1" in the H5 CRISPR Capture features. If your guides.csv still
# uses the original A/B names, ciri_step01_assignment() will fail to match
# them and return 0 assigned cells.
#
# This script reads your original guides CSV and the CellRanger
# protospacer_calls_per_cell.csv, figures out the correct name mapping,
# and writes a harmonised guides CSV you can use directly with step 01.
#
# INPUTS
# ------
#   --guides       Original guides CSV (no header: feature, type, fixed)
#   --protospacer  protospacer_calls_per_cell.csv from CellRanger output
#   --out          Output path for harmonised guides CSV
#                  Default: guides_harmonised.csv in the same dir as --guides
#
# OUTPUT
# ------
#   guides_harmonised.csv   Same format as guides CSV (no header: feature,
#                           type, fixed) but with feature names matching the
#                           H5 exactly. One row per unique CellRanger name.
#   name_mapping.csv        Full mapping table showing how each original guide
#                           name was resolved (useful for QC).
#
# USAGE
# -----
#   Rscript harmonise_guide_names.R \
#     --guides      /path/to/guides_2.csv \
#     --protospacer /path/to/protospacer_calls_per_cell.csv \
#     --out         /path/to/guides_harmonised.csv
#
# THEN use guides_harmonised.csv in step 01:
#   ciri_step01_assignment(
#     data_dir = "/path/to/data",
#     matrix   = "filtered_feature_bc_matrix.h5",
#     guides   = "guides_harmonised.csv",
#     ...
#   )
# =============================================================================

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) && idx + 1 <= length(args)) return(args[idx + 1])
  if (!is.null(default)) return(default)
  stop(paste("Missing required argument:", flag), call. = FALSE)
}

if ("--help" %in% args || "-h" %in% args) {
  cat("\nUsage: Rscript harmonise_guide_names.R",
      "--guides <file> --protospacer <file> [--out <file>]\n\n")
  quit(status = 0)
}

guides_path     <- get_arg("--guides")
protospacer_path <- get_arg("--protospacer")
out_path        <- get_arg("--out", default = NULL)

# Default output: guides_harmonised.csv next to the guides file
if (is.null(out_path)) {
  out_path <- file.path(dirname(guides_path), "guides_harmonised.csv")
}
mapping_path <- file.path(dirname(out_path), "name_mapping.csv")

cat("=======================================================================\n")
cat("  harmonise_guide_names.R\n")
cat("=======================================================================\n")
cat("  Guides CSV   :", guides_path, "\n")
cat("  Protospacer  :", protospacer_path, "\n")
cat("  Output       :", out_path, "\n")
cat("  Name mapping :", mapping_path, "\n")
cat("=======================================================================\n\n")

# ---------------------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------------------
if (!file.exists(guides_path))
  stop("guides file not found: ", guides_path, call. = FALSE)
if (!file.exists(protospacer_path))
  stop("protospacer file not found: ", protospacer_path, call. = FALSE)

guides <- read.csv(guides_path, header = FALSE, stringsAsFactors = FALSE)
if (ncol(guides) < 3)
  stop("guides CSV must have at least 3 columns: feature, type, fixed", call. = FALSE)
names(guides)[1:3] <- c("feature", "type", "fixed")

proto <- read.csv(protospacer_path, header = TRUE, stringsAsFactors = FALSE)
if (!"feature_call" %in% names(proto))
  stop("protospacer CSV must have a 'feature_call' column", call. = FALSE)

# Extract all unique guide names from the protospacer file
proto_names <- unique(unlist(strsplit(proto$feature_call, "\\|")))
proto_names <- proto_names[nchar(proto_names) > 0]
cat(sprintf("[INFO]  %d unique guide names found in protospacer CSV\n",
            length(proto_names)))
cat(sprintf("[INFO]  %d guide entries in guides CSV\n", nrow(guides)))

# ---------------------------------------------------------------------------
# Name matching
# Three rules applied in order; first match wins:
#   1. Strip trailing A or B after _number  e.g. ATF7IP_1A -> ATF7IP_1
#   2. Gene prefix (before first underscore) e.g. NTCa_1A  -> NTCa
#   3. Exact match                           e.g. MYOD_1    -> MYOD_1
# ---------------------------------------------------------------------------
find_match <- function(name, data_names) {
  # Rule 1: strip trailing A/B after _number
  m <- regmatches(name, regexpr("^.+_\\d+(?=[AB]$)", name, perl = TRUE))
  if (length(m) > 0 && m %in% data_names) return(m)

  # Rule 2: gene prefix (everything before first underscore)
  gene <- strsplit(name, "_")[[1]][1]
  if (gene %in% data_names) return(gene)

  # Rule 3: exact
  if (name %in% data_names) return(name)

  return(NA_character_)
}

guides$cellranger_name <- vapply(guides$feature, find_match,
                                  character(1), data_names = proto_names)

# ---------------------------------------------------------------------------
# Report matching results
# ---------------------------------------------------------------------------
matched   <- guides[!is.na(guides$cellranger_name), ]
unmatched <- guides[is.na(guides$cellranger_name), ]

cat(sprintf("\n[INFO]  Matched   : %d / %d guides\n", nrow(matched), nrow(guides)))

if (nrow(unmatched) > 0) {
  cat(sprintf("[WARN]  Unmatched : %d guide(s) — not present in protospacer CSV:\n",
              nrow(unmatched)))
  for (name in unmatched$feature) {
    cat(sprintf("          %s\n", name))
  }
  cat("\n  These guides will be excluded from the harmonised CSV.\n")
  cat("  Possible reasons:\n")
  cat("    - The guide was not detected in this experiment\n")
  cat("    - The guide name format is unusual — check name_mapping.csv\n\n")
}

# ---------------------------------------------------------------------------
# Build harmonised guides CSV
# One row per unique CellRanger name (deduplicate A/B replicates).
# When two original rows map to the same CellRanger name (e.g. ATF7IP_1A and
# ATF7IP_1B both -> ATF7IP_1), they must agree on type and fixed. Warn if not.
# ---------------------------------------------------------------------------
harmonised <- list()

for (cr_name in unique(matched$cellranger_name)) {
  rows <- matched[matched$cellranger_name == cr_name, ]

  # Check consistency
  if (length(unique(rows$type)) > 1)
    warning(sprintf("Guide '%s': conflicting type values (%s) — using first",
                    cr_name, paste(unique(rows$type), collapse = "/")))
  if (length(unique(rows$fixed)) > 1)
    warning(sprintf("Guide '%s': conflicting fixed values (%s) — using first",
                    cr_name, paste(unique(rows$fixed), collapse = "/")))

  harmonised[[cr_name]] <- data.frame(
    feature = cr_name,
    type    = rows$type[1],
    fixed   = rows$fixed[1],
    stringsAsFactors = FALSE
  )
}

harmonised_df <- do.call(rbind, harmonised)
rownames(harmonised_df) <- NULL

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
write.csv(harmonised_df, out_path, row.names = FALSE, quote = FALSE)
cat(sprintf("[INFO]  Harmonised guides CSV saved: %s\n", out_path))
cat(sprintf("        %d rows (from %d original, deduplicated A/B replicates)\n\n",
            nrow(harmonised_df), nrow(guides)))

# Full mapping table for reference
mapping_df <- data.frame(
  original_name    = guides$feature,
  cellranger_name  = guides$cellranger_name,
  type             = guides$type,
  fixed            = guides$fixed,
  status           = ifelse(is.na(guides$cellranger_name), "unmatched", "matched"),
  stringsAsFactors = FALSE
)
write.csv(mapping_df, mapping_path, row.names = FALSE, quote = FALSE)
cat(sprintf("[INFO]  Name mapping table saved:     %s\n\n", mapping_path))

# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
cat("=== Harmonised guide annotation ===\n")
cat(sprintf("  %-18s  %-6s  %-6s\n", "Feature", "Type", "Fixed"))
cat(sprintf("  %s\n", strrep("-", 34)))
for (i in seq_len(nrow(harmonised_df))) {
  cat(sprintf("  %-18s  %-6s  %-6s\n",
              harmonised_df$feature[i],
              harmonised_df$type[i],
              harmonised_df$fixed[i]))
}

cat("\n=== Next step ===\n")
cat(sprintf("  Use '%s' in ciri_step01_assignment():\n\n",
            basename(out_path)))
cat(sprintf("  ciri_step01_assignment(\n"))
cat(sprintf("    data_dir = \"/path/to/data\",\n"))
cat(sprintf("    matrix   = \"filtered_feature_bc_matrix.h5\",\n"))
cat(sprintf("    guides   = \"%s\",\n", basename(out_path)))
cat(sprintf("    sample   = \"your_sample\",\n"))
cat(sprintf("    strategy = 2\n"))
cat(sprintf("  )\n\n"))
