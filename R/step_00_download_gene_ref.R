#' Download Ensembl protein-coding gene reference
#'
#' Run once before the pipeline. Downloads human protein-coding genes from
#' Ensembl via biomaRt and saves a local CSV used by \code{\link{ciri_step02_filter}}.
#'
#' @param output_root Top-level output directory. Default \code{"Output"}.
#' @param sample      Experiment name used in the output folder name.
#'
#' @return Invisibly returns the path to the saved CSV.
#'
#' @section Output structure:
#' \preformatted{
#' Output/<YYMMDD>_step00_geneRef_<sample>/
#'   csv/        ensembl_protein_coding_genes.csv
#'   to_scratch/ ensembl_protein_coding_genes.csv
#' }
#'
#' @examples
#' \dontrun{
#' ciri_step00_download_ref(output_root = "Output", sample = "AB011")
#' }
#' @export
ciri_step00_download_ref <- function(output_root = "Output",
                                     sample       = "CIRI") {
  .check_bioc_pkgs("biomaRt")

  out <- make_out_dirs(output_root, "step00_geneRef", sample)

  step_banner("00", "Download Ensembl Gene Reference",
    inputs  = "Ensembl REST API (internet required)",
    outputs = file.path(out$csv, "ensembl_protein_coding_genes.csv")
  )

  options(timeout = 600)
  mirrors <- c("https://www.ensembl.org", "https://useast.ensembl.org",
               "https://uswest.ensembl.org", "https://asia.ensembl.org")
  mart <- NULL
  for (host in mirrors) {
    log_info("Trying: ", host)
    tryCatch({
      mart <- biomaRt::useMart("ENSEMBL_MART_ENSEMBL",
                               dataset = "hsapiens_gene_ensembl", host = host)
      log_info("Connected."); break
    }, error = function(e) log_warn("  Failed: ", conditionMessage(e)))
  }
  if (is.null(mart)) stop("Could not connect to any Ensembl mirror.", call. = FALSE)

  log_info("Querying protein-coding genes ...")
  ref <- biomaRt::getBM(
    mart       = mart,
    attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"),
    filters    = "biotype",
    values     = "protein_coding"
  )
  ref$mix     <- ifelse(is.na(ref$hgnc_symbol) | ref$hgnc_symbol == "",
                        ref$ensembl_gene_id, ref$hgnc_symbol)
  ref$is_mito <- ref$chromosome_name == "MT"

  out_csv <- file.path(out$csv, "ensembl_protein_coding_genes.csv")
  utils::write.csv(ref, out_csv, row.names = FALSE)
  file.copy(out_csv,
            file.path(out$to_scratch, "ensembl_protein_coding_genes.csv"),
            overwrite = TRUE)

  log_info("Saved ", nrow(ref), " genes (", sum(ref$is_mito), " mito).")
  log_info("Output folder: ", out$out_dir)
  invisible(out_csv)
}
