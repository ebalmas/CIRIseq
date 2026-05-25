#!/usr/bin/env Rscript
# =============================================================================
# CIRI Pipeline — Step 00: Download Ensembl Protein-Coding Gene Reference
# =============================================================================
# Run ONCE before the pipeline. Downloads human protein-coding genes from
# Ensembl and saves a local CSV used by Step 02.
#
# OUTPUT STRUCTURE:
#   Output/
#     <YYMMDD>_step00_geneRef_<sample>/
#       csv/    ensembl_protein_coding_genes.csv   ← passed to Step 02
#
# Usage:
#   Rscript steps/step_00_download_gene_ref.R \
#     --output_root Output \
#     --sample      AB011
#   Rscript steps/step_00_download_gene_ref.R --help
# =============================================================================

STEPS_DIR <- dirname(sys.frame(1)$ofile)
source(file.path(STEPS_DIR, "utils.R"))

SPEC <- list(
  output_root = list(flag="--output_root", required=FALSE, type="character",
                     default="Output",
                     help="Top-level output directory (will be created if needed)"),
  sample      = list(flag="--sample",      required=FALSE, type="character",
                     default="CIRI",
                     help="Sample / experiment name used in the output folder name")
)

main <- function() {
  p   <- parse_args(SPEC)
  out <- make_out_dirs(p$output_root, "step00_geneRef", p$sample)

  step_banner("00", "Download Ensembl Gene Reference",
    inputs  = c("Ensembl REST API (internet required)"),
    outputs = c(file.path(out$csv, "ensembl_protein_coding_genes.csv"))
  )

  suppressMessages(library(biomaRt))
  options(timeout = 600)

  mirrors <- c("https://www.ensembl.org", "https://useast.ensembl.org",
               "https://uswest.ensembl.org", "https://asia.ensembl.org")
  mart <- NULL
  for (host in mirrors) {
    log_info("Trying: ", host)
    tryCatch({
      mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = host)
      log_info("Connected to: ", host); break
    }, error = function(e) log_warn("  Failed: ", conditionMessage(e)))
  }
  if (is.null(mart)) stop("Could not connect to any Ensembl mirror.")

  log_info("Querying protein-coding genes ...")
  ref <- getBM(mart = mart,
               attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"),
               filters    = "biotype", values = "protein_coding")
  ref$mix     <- ifelse(is.na(ref$hgnc_symbol) | ref$hgnc_symbol == "",
                        ref$ensembl_gene_id, ref$hgnc_symbol)
  ref$is_mito <- ref$chromosome_name == "MT"

  out_csv <- file.path(out$csv, "ensembl_protein_coding_genes.csv")
  write.csv(ref, out_csv, row.names = FALSE)

  # Also copy to to_scratch for Step 02 to pick up automatically
  file.copy(out_csv, file.path(out$to_scratch, "ensembl_protein_coding_genes.csv"))

  log_info("Saved ", nrow(ref), " genes (", sum(ref$is_mito), " mito) to:")
  log_info("  ", out_csv)
  log_info("Step 00 complete. to_scratch/ ready for Step 02.")
}

main()
