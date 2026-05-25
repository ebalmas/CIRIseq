# =============================================================================
# CIRI — Scratch / staging area helpers
# =============================================================================

#' List files in a step's to_scratch/ folder
#'
#' Shows what a completed step has produced and is ready to be promoted.
#'
#' @param step_tag    Step tag, e.g. \code{"step01_assignment"}.
#' @param sample      Experiment name.
#' @param output_root Top-level output directory. Default \code{"Output"}.
#'
#' @return Invisibly returns a character vector of file paths.
#' @examples
#' \dontrun{
#' ciri_list_scratch("step01_assignment", sample = "AB011")
#' }
#' @export
ciri_list_scratch <- function(step_tag,
                               sample       = "CIRI",
                               output_root  = "Output") {
  to_scratch <- resolve_scratch(NULL, output_root, step_tag, sample)
  files      <- list.files(to_scratch, full.names = TRUE, recursive = TRUE)
  if (!length(files)) {
    message("  (to_scratch/ is empty for ", step_tag, " / ", sample, ")")
  } else {
    message("Files ready in to_scratch/ for ", step_tag, " / ", sample, ":")
    for (f in files) message("  ", f)
  }
  invisible(files)
}

#' Promote to_scratch/ output to scratch/ staging area
#'
#' After inspecting a step's output and deciding to proceed, call this to
#' copy the files from \code{to_scratch/} into the project-level
#' \code{scratch/} folder so the next step can read them.
#'
#' This is the deliberate checkpoint: nothing moves to \code{scratch/}
#' automatically. You decide when you are happy with the output.
#'
#' @param step_tag    Step tag of the step you just finished,
#'                    e.g. \code{"step01_assignment"}.
#' @param sample      Experiment name.
#' @param output_root Top-level output directory. Default \code{"Output"}.
#' @param scratch_dir Project-level scratch directory. Default \code{"scratch"}.
#' @param overwrite   If \code{TRUE}, overwrite existing files in scratch/.
#'                    Default \code{FALSE} (safe).
#'
#' @return Invisibly returns the scratch directory path.
#'
#' @details
#' \strong{Workflow:}
#' \preformatted{
#'   ciri_step01_assignment(...)
#'   # → inspect Output/<date>_step01.../plots/ and stats/
#'   ciri_promote_scratch("step01_assignment", sample = "AB011")
#'   # → files now in scratch/  (Step 02 will read from here)
#'   ciri_step02_filter(...)
#' }
#'
#' @examples
#' \dontrun{
#' ciri_promote_scratch("step01_assignment", sample = "AB011")
#' }
#' @export
ciri_promote_scratch <- function(step_tag,
                                  sample       = "CIRI",
                                  output_root  = "Output",
                                  scratch_dir  = "scratch",
                                  overwrite    = FALSE) {

  to_scratch <- resolve_scratch(NULL, output_root, step_tag, sample)
  files      <- list.files(to_scratch, full.names = TRUE, recursive = TRUE)

  if (!length(files)) stop("Nothing in to_scratch/ for ", step_tag, " / ", sample)

  dir.create(scratch_dir, showWarnings = FALSE, recursive = TRUE)

  message(sprintf("Promoting %d file(s) from:\n  %s\nto:\n  %s\n",
                  length(files), to_scratch, scratch_dir))

  for (f in files) {
    dest <- file.path(scratch_dir, basename(f))
    if (file.exists(dest) && !overwrite) {
      message("  [SKIP — already exists] ", basename(f),
              "  (use overwrite = TRUE to replace)")
      next
    }
    file.copy(f, dest, overwrite = overwrite)
    message("  [OK] ", basename(f))
  }

  message("\nscratch/ is ready. Proceed with the next step.")
  invisible(scratch_dir)
}

#' Show the current contents of scratch/
#'
#' @param scratch_dir Path to the scratch directory. Default \code{"scratch"}.
#' @return Invisibly returns a character vector of file paths.
#' @examples
#' \dontrun{
#' ciri_scratch_status()
#' }
#' @export
ciri_scratch_status <- function(scratch_dir = "scratch") {
  if (!dir.exists(scratch_dir)) {
    message("scratch/ does not exist yet. Run ciri_promote_scratch() first.")
    return(invisible(character(0)))
  }
  files <- list.files(scratch_dir, full.names = TRUE, recursive = TRUE)
  if (!length(files)) {
    message("scratch/ is empty.")
  } else {
    message("Current scratch/ contents:")
    for (f in files) message("  ", basename(f), "  (", file.size(f), " bytes)")
  }
  invisible(files)
}

#' Clear the scratch/ staging area
#'
#' Use before starting a new step to avoid mixing outputs from different runs.
#'
#' @param scratch_dir Path to the scratch directory. Default \code{"scratch"}.
#' @param confirm     Must be set to \code{TRUE} to actually delete. Safety guard.
#' @return Invisibly \code{NULL}.
#' @examples
#' \dontrun{
#' ciri_clear_scratch(confirm = TRUE)
#' }
#' @export
ciri_clear_scratch <- function(scratch_dir = "scratch", confirm = FALSE) {
  if (!confirm)
    stop("Set confirm = TRUE to clear scratch/. This deletes all files in that folder.")
  files <- list.files(scratch_dir, full.names = TRUE)
  file.remove(files)
  message("scratch/ cleared (", length(files), " files removed).")
  invisible(NULL)
}
