# =============================================================================
# CIRI — Dependency installer helper
# =============================================================================

#' Install all optional dependencies for the CIRI pipeline
#'
#' CIRI itself only requires lightweight CRAN packages and installs quickly.
#' The heavy bioinformatics packages (\code{Seurat}, \code{monocle3},
#' \code{hdf5r}, \code{biomaRt}) are optional and only loaded when a specific
#' step needs them. This function installs all of them at once when you are
#' ready to run the full pipeline.
#'
#' You can also install subsets using the \code{steps} argument if you only
#' plan to run part of the pipeline.
#'
#' @param steps Character vector of steps to install dependencies for.
#'   Use \code{"all"} (default) to install everything.
#'   Options: \code{"00"}, \code{"01"}, \code{"02"}, \code{"03-08"}, \code{"all"}.
#' @param upgrade Passed to \code{BiocManager::install()}. Default \code{"never"}
#'   (do not upgrade already-installed packages).
#'
#' @return Invisibly \code{NULL}. Called for side effects.
#'
#' @section What gets installed:
#' \describe{
#'   \item{Step 00}{biomaRt (Bioconductor)}
#'   \item{Step 01}{hdf5r (CRAN) — also needs the HDF5 system library, see README}
#'   \item{Step 02}{hdf5r, Seurat, data.table (CRAN)}
#'   \item{Steps 03-08}{monocle3 (GitHub), igraph (CRAN)}
#' }
#'
#' @examples
#' \dontrun{
#' # Install everything
#' CIRI::install_dependencies()
#'
#' # Install only what you need for Steps 01-02
#' CIRI::install_dependencies(steps = c("01", "02"))
#'
#' # Install only Monocle3 and trajectory dependencies
#' CIRI::install_dependencies(steps = "03-08")
#' }
#'
#' @export
install_dependencies <- function(steps = "all", upgrade = "never") {

  # Map steps to required packages
  step_pkgs <- list(
    "00"    = list(cran = character(0),
                   bioc = "biomaRt",
                   gh   = character(0)),
    "01"    = list(cran = "hdf5r",
                   bioc = character(0),
                   gh   = character(0)),
    "02"    = list(cran = c("hdf5r", "Seurat", "data.table"),
                   bioc = character(0),
                   gh   = character(0)),
    "03-08" = list(cran = "igraph",
                   bioc = character(0),
                   gh   = "cole-trapnell-lab/monocle3")
  )

  if (identical(steps, "all")) steps <- names(step_pkgs)

  # Collect unique packages across requested steps
  to_cran <- unique(unlist(lapply(step_pkgs[steps], `[[`, "cran")))
  to_bioc <- unique(unlist(lapply(step_pkgs[steps], `[[`, "bioc")))
  to_gh   <- unique(unlist(lapply(step_pkgs[steps], `[[`, "gh")))

  # Filter out already-installed packages
  to_cran <- to_cran[!vapply(to_cran, requireNamespace, logical(1), quietly = TRUE)]
  to_bioc <- to_bioc[!vapply(to_bioc, requireNamespace, logical(1), quietly = TRUE)]
  to_gh_names <- gsub(".*/", "", to_gh)  # e.g. "monocle3" from "cole-trapnell-lab/monocle3"
  to_gh   <- to_gh[!vapply(to_gh_names,  requireNamespace, logical(1), quietly = TRUE)]

  if (length(to_cran) == 0 && length(to_bioc) == 0 && length(to_gh) == 0) {
    message("All requested dependencies are already installed.")
    return(invisible(NULL))
  }

  # ---- CRAN -----------------------------------------------------------------
  if (length(to_cran) > 0) {
    message("\nInstalling from CRAN: ", paste(to_cran, collapse = ", "))
    if ("hdf5r" %in% to_cran) {
      message(
        "\n  NOTE: hdf5r requires the HDF5 system library.\n",
        "  macOS:         brew install hdf5\n",
        "  Ubuntu/Debian: sudo apt-get install libhdf5-dev\n",
        "  Fedora/RHEL:   sudo dnf install hdf5-devel\n",
        "  Install that first if hdf5r fails.\n"
      )
    }
    utils::install.packages(to_cran)
  }

  # ---- Bioconductor ---------------------------------------------------------
  if (length(to_bioc) > 0) {
    message("\nInstalling from Bioconductor: ", paste(to_bioc, collapse = ", "))
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      message("  Installing BiocManager first ...")
      utils::install.packages("BiocManager")
    }
    BiocManager::install(to_bioc, upgrade = upgrade, ask = FALSE)
  }

  # ---- GitHub ---------------------------------------------------------------
  if (length(to_gh) > 0) {
    message("\nInstalling from GitHub: ", paste(to_gh, collapse = ", "))
    if (!requireNamespace("remotes", quietly = TRUE)) {
      message("  Installing remotes first ...")
      utils::install.packages("remotes")
    }
    for (repo in to_gh) {
      message("  Installing: ", repo)
      remotes::install_github(repo, upgrade = upgrade)
    }
  }

  message("\nDone. Restart R if prompted.")
  invisible(NULL)
}

#' Check which optional dependencies are installed
#'
#' Prints a summary of which pipeline steps are ready to run based on
#' currently installed packages.
#'
#' @return Invisibly returns a named logical vector.
#' @examples
#' \dontrun{
#' CIRI::check_dependencies()
#' }
#' @export
check_dependencies <- function() {

  pkgs <- list(
    "Step 00 — Ensembl reference"  = "biomaRt",
    "Step 01 — Guide assignment"   = "hdf5r",
    "Step 02 — Filter"             = c("hdf5r", "Seurat", "data.table"),
    "Steps 03-08 — Monocle3"       = c("monocle3", "igraph")
  )

  message("\nCIRI dependency status:")
  message(strrep("-", 50))

  status <- vapply(names(pkgs), function(step) {
    needed  <- pkgs[[step]]
    missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
    ok      <- length(missing) == 0
    icon    <- if (ok) "[OK]  " else "[MISS]"
    msg     <- if (ok) "ready" else paste("missing:", paste(missing, collapse = ", "))
    message(sprintf("  %s  %-36s %s", icon, step, msg))
    ok
  }, logical(1))

  message(strrep("-", 50))
  n_ready <- sum(status)
  if (n_ready == length(status)) {
    message("  All dependencies installed. Pipeline is ready to run.\n")
  } else {
    message(sprintf("  %d/%d step groups ready. Run install_dependencies() to install missing packages.\n",
                    n_ready, length(status)))
  }

  invisible(status)
}
