# =============================================================================
# CIRI Pipeline — Shared Utilities
# Source this at the top of every step: source(file.path(STEPS_DIR, "utils.R"))
# =============================================================================

# ---------------------------------------------------------------------------
# Output folder structure (mirrors catcheR convention)
#
#   Output/
#     <YYMMDD>_<step>_<sample>/
#       csv/          tables
#       plots/        PDF figures
#       stats/        text summaries / logs
#       R_objects/    .RData / .rds objects
#       to_scratch/   files needed by the NEXT step
#
# Each step returns a list: $out_dir (the dated root) and sub-paths.
# ---------------------------------------------------------------------------

#' Create the dated output folder tree for one step.
#'
#' @param output_root  Top-level Output/ directory (absolute or relative).
#' @param step_tag     Short step identifier, e.g. "step01_assignment".
#' @param sample       Sample / experiment name, e.g. "AB011".
#' @return Named list: out_dir, csv, plots, stats, R_objects, to_scratch
make_out_dirs <- function(output_root, step_tag, sample = "CIRI") {
  date_str <- format(Sys.Date(), "%y%m%d")
  folder   <- file.path(output_root, paste(date_str, step_tag, sample, sep = "_"))

  dirs <- list(
    out_dir   = folder,
    csv       = file.path(folder, "csv"),
    plots     = file.path(folder, "plots"),
    stats     = file.path(folder, "stats"),
    R_objects = file.path(folder, "R_objects"),
    to_scratch = file.path(folder, "to_scratch")
  )
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dirs
}

# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------
# spec <- list(
#   dir = list(flag="--dir", required=TRUE, type="character", help="Data dir")
# )
# p <- parse_args(spec)

parse_args <- function(spec) {
  raw <- commandArgs(trailingOnly = TRUE)

  if ("--help" %in% raw || "-h" %in% raw) {
    cat("\nArguments:\n")
    for (nm in names(spec)) {
      s       <- spec[[nm]]
      req_tag <- if (isTRUE(s$required)) " [required]" else ""
      def_tag <- if (!is.null(s$default)) paste0(" (default: ", s$default, ")") else ""
      cat(sprintf("  %-22s %s%s%s\n", s$flag, s$help, req_tag, def_tag))
    }
    cat("\n")
    quit(status = 0)
  }

  result <- list()
  for (nm in names(spec)) {
    s <- spec[[nm]]
    if (!is.null(s$default)) result[[nm]] <- s$default
  }

  i <- 1
  while (i <= length(raw)) {
    matched <- FALSE
    for (nm in names(spec)) {
      if (raw[i] == spec[[nm]]$flag) {
        val_raw    <- raw[i + 1]
        result[[nm]] <- switch(spec[[nm]]$type,
          "integer" = as.integer(val_raw),
          "numeric" = as.numeric(val_raw),
          "logical" = toupper(val_raw) %in% c("TRUE", "T", "YES", "1"),
          val_raw
        )
        i       <- i + 2
        matched <- TRUE
        break
      }
    }
    if (!matched) i <- i + 1
  }

  missing_flags <- names(spec)[sapply(names(spec), function(nm)
    isTRUE(spec[[nm]]$required) && is.null(result[[nm]]))]
  if (length(missing_flags) > 0)
    stop(paste("Missing required arguments:",
               paste(sapply(missing_flags, function(nm) spec[[nm]]$flag), collapse = ", "),
               "\nRun with --help for usage."))

  result
}

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
log_info  <- function(...) message("[INFO]  ", ...)
log_warn  <- function(...) message("[WARN]  ", ...)
log_error <- function(...) message("[ERROR] ", ...)

# ---------------------------------------------------------------------------
# Step banner — printed at runtime so you always see inputs/outputs
# ---------------------------------------------------------------------------
step_banner <- function(n, title, inputs = NULL, outputs = NULL) {
  sep <- strrep("=", 72)
  message(sep)
  message(sprintf("  CIRI Pipeline — Step %s: %s", n, title))
  message(sprintf("  %s", Sys.time()))
  message(sep)
  if (!is.null(inputs)) {
    message("  INPUTS :")
    for (x in inputs) message("    <  ", x)
  }
  if (!is.null(outputs)) {
    message("  OUTPUTS:")
    for (x in outputs) message("    >  ", x)
  }
  message(sep)
}

# ---------------------------------------------------------------------------
# Load .RData returning the first object (any variable name)
# ---------------------------------------------------------------------------
load_rdata <- function(path) {
  if (!file.exists(path)) stop(paste("File not found:", path))
  env <- new.env()
  load(path, envir = env)
  env[[ls(env)[1]]]
}

# ---------------------------------------------------------------------------
# Assert a file exists — clear error + hint
# ---------------------------------------------------------------------------
assert_file <- function(path, hint = NULL) {
  if (!file.exists(path)) {
    msg <- paste("Required file not found:", path)
    if (!is.null(hint)) msg <- paste0(msg, "\nHint: ", hint)
    stop(msg)
  }
  invisible(path)
}

# ---------------------------------------------------------------------------
# Resolve to_scratch path from a previous step's output root.
# Looks for Output/<date>_<step_tag>_<sample>/to_scratch/
# If multiple dated folders exist, takes the most recent.
# ---------------------------------------------------------------------------
find_scratch <- function(output_root, step_tag, sample = "CIRI") {
  pattern <- paste0("*_", step_tag, "_", sample)
  hits    <- list.dirs(output_root, recursive = FALSE)
  hits    <- hits[grepl(paste0("_", step_tag, "_", sample, "$"), hits)]
  if (length(hits) == 0)
    stop(paste0("No output folder found for step '", step_tag,
                "' sample '", sample, "' under ", output_root,
                "\nRun the previous step first."))
  scratch <- file.path(sort(hits, decreasing = TRUE)[1], "to_scratch")
  if (!dir.exists(scratch))
    stop(paste("to_scratch/ folder not found in:", sort(hits)[1]))
  scratch
}
