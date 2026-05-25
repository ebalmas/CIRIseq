#!/usr/bin/env Rscript
# =============================================================================
# CIRI Pipeline - Step 00: Download Ensembl Protein-Coding Gene Reference
# =============================================================================
# Run this ONCE before the pipeline. It saves a local CSV that anno_filter.R
# will pick up automatically, avoiding repeated Ensembl queries.
#
# Usage:
#   Rscript 00_download_gene_ref.R [output_path]
#
# Arguments:
#   output_path  Where to save the CSV.
#                Default: ensembl_protein_coding_genes.csv (working dir)
# =============================================================================

suppressMessages({
  library(biomaRt)
})

args <- commandArgs(trailingOnly = TRUE)
out_file <- if (length(args) >= 1) args[1] else "ensembl_protein_coding_genes.csv"

message("=== CIRI - Ensembl Reference Downloader ===")
message(paste("Output file:", out_file))

# Try a series of Ensembl mirrors in order of reliability
hosts <- c(
  "https://www.ensembl.org",
  "https://useast.ensembl.org",
  "https://uswest.ensembl.org",
  "https://asia.ensembl.org"
)

mart <- NULL
for (host in hosts) {
  message(paste("Trying host:", host))
  tryCatch({
    mart <- useMart(
      biomart = "ENSEMBL_MART_ENSEMBL",
      dataset = "hsapiens_gene_ensembl",
      host    = host
    )
    message(paste("Connected to:", host))
    break
  }, error = function(e) {
    message(paste("  Failed:", conditionMessage(e)))
  })
}

if (is.null(mart)) {
  stop("Could not connect to any Ensembl mirror. Check your internet connection.")
}

message("Downloading protein-coding gene list ...")
options(timeout = 600)

all_coding_genes <- getBM(
  mart       = mart,
  attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name"),
  filters    = "biotype",
  values     = "protein_coding"
)

# Fill missing HGNC symbols with Ensembl IDs
all_coding_genes$mix <- ifelse(
  is.na(all_coding_genes$hgnc_symbol) | all_coding_genes$hgnc_symbol == "",
  all_coding_genes$ensembl_gene_id,
  all_coding_genes$hgnc_symbol
)

# Flag mitochondrial genes (chrMT) — useful for downstream reference
all_coding_genes$is_mitochondrial <- all_coding_genes$chromosome_name == "MT"

write.csv(all_coding_genes, out_file, row.names = FALSE)

n_total <- nrow(all_coding_genes)
n_mito  <- sum(all_coding_genes$is_mitochondrial, na.rm = TRUE)
message(paste("Done. Downloaded", n_total, "protein-coding genes (", n_mito, "mitochondrial)."))
message(paste("Saved to:", out_file))
