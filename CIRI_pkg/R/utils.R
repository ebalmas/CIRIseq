# =============================================================================
# CIRI — Internal utilities (not exported)
# =============================================================================

# ---------------------------------------------------------------------------
# Output folder structure
# ---------------------------------------------------------------------------

#' Create dated output subdirectory tree for one step.
#'
#' Structure:
#'   <output_root>/<YYMMDD>_<step_tag>_<sample>/
#'     csv/  plots/  stats/  R_objects/  to_scratch/
#'
#' @param output_root Top-level Output/ directory.
#' @param step_tag    Short step id, e.g. "step01_assignment".
#' @param sample      Experiment name, e.g. "AB011".
#' @return Named list of directory paths.
#' @keywords internal
make_out_dirs <- function(output_root, step_tag, sample = "CIRI") {
  date_str <- format(Sys.Date(), "%y%m%d")
  folder   <- file.path(output_root,
                        paste(date_str, step_tag, sample, sep = "_"))
  dirs <- list(
    out_dir    = folder,
    csv        = file.path(folder, "csv"),
    plots      = file.path(folder, "plots"),
    stats      = file.path(folder, "stats"),
    R_objects  = file.path(folder, "R_objects"),
    to_scratch = file.path(folder, "to_scratch")
  )
  for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  dirs
}

# ---------------------------------------------------------------------------
# scratch/ resolution
# ---------------------------------------------------------------------------

#' Find the scratch/ folder for a given step and sample.
#'
#' Reads from <output_root>/<date>_<step_tag>_<sample>/to_scratch/
#' unless the user overrides with an explicit scratch path.
#'
#' @param scratch     Explicit path (if not NULL, returned as-is).
#' @param output_root Top-level Output/ directory.
#' @param step_tag    Step tag of the PREVIOUS step.
#' @param sample      Sample name.
#' @return Absolute path to the to_scratch/ directory.
#' @keywords internal
resolve_scratch <- function(scratch = NULL, output_root, step_tag, sample) {
  if (!is.null(scratch)) {
    if (!dir.exists(scratch))
      stop("scratch directory not found: ", scratch)
    return(scratch)
  }
  hits    <- list.dirs(output_root, recursive = FALSE)
  pattern <- paste0("_", step_tag, "_", sample, "$")
  hits    <- hits[grepl(pattern, hits)]
  if (!length(hits))
    stop("No output folder found for step '", step_tag, "' sample '", sample,
         "' under ", output_root,
         "\nRun the previous step first, or pass scratch= explicitly.")
  to_scratch <- file.path(sort(hits, decreasing = TRUE)[1], "to_scratch")
  if (!dir.exists(to_scratch))
    stop("to_scratch/ not found in: ", sort(hits, decreasing = TRUE)[1])
  to_scratch
}

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
log_info  <- function(...) message("[INFO]  ", ...)
log_warn  <- function(...) message("[WARN]  ", ...)

# ---------------------------------------------------------------------------
# Step banner
# ---------------------------------------------------------------------------
step_banner <- function(n, title, inputs = NULL, outputs = NULL) {
  sep <- strrep("=", 72)
  message(sep)
  message(sprintf("  CIRI — Step %s: %s", n, title))
  message(sprintf("  %s", Sys.time()))
  message(sep)
  if (!is.null(inputs))  { message("  INPUTS :"); for (x in inputs)  message("    <  ", x) }
  if (!is.null(outputs)) { message("  OUTPUTS:"); for (x in outputs) message("    >  ", x) }
  message(sep)
}

# ---------------------------------------------------------------------------
# Load .RData — returns first object regardless of variable name
# ---------------------------------------------------------------------------
load_rdata <- function(path) {
  if (!file.exists(path)) stop("File not found: ", path)
  env <- new.env()
  load(path, envir = env)
  env[[ls(env)[1]]]
}

# ---------------------------------------------------------------------------
# Assert file exists
# ---------------------------------------------------------------------------
assert_file <- function(path, hint = NULL) {
  if (!file.exists(path)) {
    msg <- paste("Required file not found:", path)
    if (!is.null(hint)) msg <- paste0(msg, "\nHint: ", hint)
    stop(msg)
  }
  invisible(path)
}
