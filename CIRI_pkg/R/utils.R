# =============================================================================
# CIRI — Internal utilities (not exported)
# =============================================================================

# ---------------------------------------------------------------------------
# requireNamespace guard — called at the top of each step function
# to give a clear error if a heavy dependency is not installed
# ---------------------------------------------------------------------------
.check_pkg <- function(pkg, install_hint) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Package '", pkg, "' is required for this step but is not installed.\n",
      "Install it with: ", install_hint,
      call. = FALSE
    )
  }
}

.check_bioc_pkgs <- function(...) {
  pkgs <- c(...)
  for (pkg in pkgs) {
    .check_pkg(pkg, paste0("BiocManager::install('", pkg, "')"))
  }
}

.check_cran_pkgs <- function(...) {
  pkgs <- c(...)
  for (pkg in pkgs) {
    .check_pkg(pkg, paste0("install.packages('", pkg, "')"))
  }
}

# ---------------------------------------------------------------------------
# Output folder structure
# ---------------------------------------------------------------------------

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

#' @keywords internal
resolve_scratch <- function(scratch = NULL, output_root, step_tag, sample) {
  if (!is.null(scratch)) {
    if (!dir.exists(scratch))
      stop("scratch directory not found: ", scratch, call. = FALSE)
    return(scratch)
  }
  hits    <- list.dirs(output_root, recursive = FALSE)
  pattern <- paste0("_", step_tag, "_", sample, "$")
  hits    <- hits[grepl(pattern, hits)]
  if (!length(hits))
    stop(
      "No output folder found for step '", step_tag, "' sample '", sample,
      "' under ", output_root,
      "\nRun the previous step first, or pass scratch= explicitly.",
      call. = FALSE
    )
  to_scratch <- file.path(sort(hits, decreasing = TRUE)[1], "to_scratch")
  if (!dir.exists(to_scratch))
    stop("to_scratch/ not found in: ", sort(hits, decreasing = TRUE)[1],
         call. = FALSE)
  to_scratch
}

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
log_info <- function(...) message("[INFO]  ", ...)
log_warn <- function(...) message("[WARN]  ", ...)

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
  if (!file.exists(path)) stop("File not found: ", path, call. = FALSE)
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
    stop(msg, call. = FALSE)
  }
  invisible(path)
}
