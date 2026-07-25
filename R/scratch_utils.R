# =============================================================================
# CIRI — Scratch / staging area helpers
# =============================================================================

#' List files in a step's to_scratch/ folder
#'
#' @param step_tag    Step tag, e.g. \code{"step01_assignment"}.
#' @param sample      Experiment name.
#' @param output_root Top-level output directory. Default \code{"Output"}.
#' @return Invisibly returns a character vector of file paths.
#' @examples
#' \dontrun{
#' ciri_list_scratch("step01_assignment", sample = "AB011")
#' }
#' @export
ciri_list_scratch <- function(step_tag, sample = "CIRI", output_root = "Output") {
  to_scratch <- resolve_scratch(NULL, output_root, step_tag, sample)
  files      <- list.files(to_scratch, full.names = TRUE, recursive = TRUE)
  if (!length(files)) message("  (to_scratch/ is empty for ", step_tag, " / ", sample, ")")
  else { message("Files ready in to_scratch/ for ", step_tag, " / ", sample, ":"); for (f in files) message("  ", f) }
  invisible(files)
}

#' Promote to_scratch/ output to scratch/ staging area
#'
#' After inspecting a step's output and deciding to proceed, call this to
#' copy the files from \code{to_scratch/} into the project-level
#' \code{scratch/} folder so the next step can read them.
#'
#' @param step_tag    Step tag, e.g. \code{"step01_assignment"}.
#' @param sample      Experiment name.
#' @param output_root Top-level output directory. Default \code{"Output"}.
#' @param scratch_dir Project-level scratch directory. Default \code{"scratch"}.
#' @param overwrite   Overwrite existing files. Default \code{FALSE}.
#' @return Invisibly returns the scratch directory path.
#' @export
ciri_promote_scratch <- function(step_tag, sample = "CIRI",
                                  output_root = "Output", scratch_dir = "scratch",
                                  overwrite = FALSE) {
  to_scratch <- resolve_scratch(NULL, output_root, step_tag, sample)
  files      <- list.files(to_scratch, full.names = TRUE, recursive = TRUE)
  if (!length(files)) stop("Nothing in to_scratch/ for ", step_tag, " / ", sample)
  dir.create(scratch_dir, showWarnings = FALSE, recursive = TRUE)
  message(sprintf("Promoting %d file(s) from:\n  %s\nto:\n  %s\n", length(files), to_scratch, scratch_dir))
  for (f in files) {
    dest <- file.path(scratch_dir, basename(f))
    if (file.exists(dest) && !overwrite) { message("  [SKIP] ", basename(f)); next }
    file.copy(f, dest, overwrite = overwrite)
    message("  [OK]   ", basename(f))
  }
  message("\nscratch/ is ready.")
  invisible(scratch_dir)
}

#' Show the current contents of scratch/
#'
#' @param scratch_dir Path to the scratch directory. Default \code{"scratch"}.
#' @return Invisibly returns a character vector of file paths.
#' @export
ciri_scratch_status <- function(scratch_dir = "scratch") {
  if (!dir.exists(scratch_dir)) { message("scratch/ does not exist yet."); return(invisible(character(0))) }
  files <- list.files(scratch_dir, full.names = TRUE, recursive = TRUE)
  if (!length(files)) message("scratch/ is empty.")
  else { message("Current scratch/ contents:"); for (f in files) message("  ", basename(f), "  (", file.size(f), " bytes)") }
  invisible(files)
}

#' Clear the scratch/ staging area
#'
#' @param scratch_dir Path to the scratch directory. Default \code{"scratch"}.
#' @param confirm     Must be \code{TRUE} to actually delete.
#' @return Invisibly \code{NULL}.
#' @export
ciri_clear_scratch <- function(scratch_dir = "scratch", confirm = FALSE) {
  if (!confirm) stop("Set confirm = TRUE to clear scratch/.")
  files <- list.files(scratch_dir, full.names = TRUE)
  file.remove(files)
  message("scratch/ cleared (", length(files), " files removed).")
  invisible(NULL)
}

#' Copy standalone QC scripts to your working directory
#'
#' Copies \code{harmonise_guide_names.R}, \code{ciri_step02_annotate.R}, and
#' \code{ciri_step02_filter.R} from the package \code{inst/scripts/} folder
#' into your working directory so you can run them directly.
#'
#' These scripts handle the Seurat QC workflow (step 02) which is kept
#' outside the main package so you can substitute your own QC pipeline.
#'
#' @param dest      Destination directory. Default: current working directory.
#' @param overwrite Overwrite existing files. Default \code{FALSE}.
#' @return Invisibly returns paths of copied files.
#' @examples
#' \dontrun{
#' ciri_copy_scripts()
#' ciri_copy_scripts(dest = "scripts/")
#' }
#' @export
ciri_copy_scripts <- function(dest = ".", overwrite = FALSE) {
  scripts_dir <- system.file("scripts", package = "CIRI")
  if (!nchar(scripts_dir)) stop("Could not find inst/scripts/ in CIRI package.", call. = FALSE)
  scripts <- list.files(scripts_dir, pattern = "\\.R$", full.names = TRUE)
  if (!length(scripts)) stop("No scripts found in inst/scripts/.", call. = FALSE)
  dir.create(dest, showWarnings = FALSE, recursive = TRUE)
  copied <- character(0)
  for (s in scripts) {
    target <- file.path(dest, basename(s))
    if (file.exists(target) && !overwrite) { message("[SKIP] ", basename(s)); next }
    file.copy(s, target, overwrite = overwrite)
    message("[OK]   ", basename(s), " -> ", target)
    copied <- c(copied, target)
  }
  message("\nEdit INTERACTIVE_PARAMS at the top of each script, then source() or Rscript.")
  invisible(copied)
}
