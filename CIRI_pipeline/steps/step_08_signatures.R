#!/usr/bin/env Rscript
# =============================================================================
# CIRI Pipeline — Step 08: Signature Scoring
# =============================================================================
# INPUT  (auto-resolved from Step 06 to_scratch/):
#   processed_cds_<group>.RData
#
# OUTPUT STRUCTURE:
#   Output/
#     <YYMMDD>_step08_signatures_<sample>/
#       csv/    signature_values/<sig>.csv   per-cell scores
#               signature_summary.csv        mean per perturbation × signature
#       plots/  umap_<sig>.pdf
#       stats/  signature_coverage.txt       genes found / total per signature
#
# To add custom gene sets: edit the SIGNATURES list below.
#
# Usage:
#   Rscript steps/step_08_signatures.R \
#     --sample AB011 \
#     --group  muscle
#   Rscript steps/step_08_signatures.R --help
# =============================================================================

STEPS_DIR <- dirname(sys.frame(1)$ofile)
source(file.path(STEPS_DIR, "utils.R"))

SPEC <- list(
  output_root = list(flag="--output_root", required=FALSE, type="character",
                     default="Output",  help="Top-level output directory"),
  sample      = list(flag="--sample",      required=FALSE, type="character",
                     default="CIRI",    help="Sample / experiment name"),
  scratch     = list(flag="--scratch",     required=FALSE, type="character",
                     default=NULL,      help="Path to Step 06 to_scratch/ (auto-resolved if omitted)"),
  group       = list(flag="--group",       required=TRUE,  type="character",
                     help="Group name matching Step 06 --group")
)

# ── Edit gene sets here ────────────────────────────────────────────────────
SIGNATURES <- list(
  cellcycle_s = c(
    "MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2","MCM6",
    "CDCA7","DTL","PRIM1","UHRF1","HELLS","RFC2","RPA2","NASP","RAD51",
    "CHEK1","ORC6","MCM3","MCM7","BRIP1","E2F8"),
  cellcycle_g2m = c(
    "HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80","CKS2",
    "NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","SMC4","CCNB2","AURKB",
    "BUB1","KIF11","CDC20","TTK"),
  sarcomere_core = c(
    "MYH3","MYH8","MYBPC1","MYBPC2","MYL1","MYL4","TNNI1","TNNI2",
    "TNNT1","TNNT2","TNNC1","TNNC2","TPM1","TPM2","ACTA1","ACTC1","NEB","NEBL"),
  myogenic_commitment = c(
    "MYOD1","MYOG","MYF5","MYF6","PAX7","MEF2A","MEF2C","MEF2D",
    "CDH15","NCAM1","DES","VIM","NES"),
  pluripotency = c(
    "POU5F1","SOX2","NANOG","KLF4","MYC","LIN28A","DPPA4","DPPA5",
    "SALL4","UTF1","DNMT3L")
)
# ──────────────────────────────────────────────────────────────────────────

main <- function() {
  suppressMessages({
    library(monocle3); library(dplyr); library(ggplot2); library(viridis)
  })

  p   <- parse_args(SPEC)
  out <- make_out_dirs(p$output_root, "step08_signatures", p$sample)
  dir.create(file.path(out$csv, "signature_values"), showWarnings=FALSE)

  scratch06 <- if (!is.null(p$scratch)) p$scratch else
    find_scratch(p$output_root, "step06_trajectory", p$sample)
  cds_path  <- file.path(scratch06, paste0("processed_cds_", p$group, ".RData"))

  step_banner("08", "Signature Scoring",
    inputs  = c(cds_path),
    outputs = c(file.path(out$csv, "signature_values/"),
                out$plots,
                file.path(out$csv, "signature_summary.csv"),
                file.path(out$stats, "signature_coverage.txt"))
  )
  assert_file(cds_path, hint = "Run step_06_subclusters_trajectory.R first.")

  cds     <- load_rdata(cds_path)
  present <- rownames(cds)

  summary_rows <- list()
  coverage_log <- c()

  for (sig_name in names(SIGNATURES)) {
    genes <- intersect(SIGNATURES[[sig_name]], present)
    n_total <- length(SIGNATURES[[sig_name]])
    n_found <- length(genes)
    coverage_log <- c(coverage_log,
      sprintf("%-30s %d / %d genes found (%.0f%%)", sig_name, n_found, n_total, n_found/n_total*100))
    log_info(tail(coverage_log, 1))

    if (n_found < 3) { log_warn("Skipping '", sig_name, "' — fewer than 3 genes found."); next }

    norm_sub <- normalized_counts(cds)[genes, , drop=FALSE]
    z_mat    <- t(scale(t(as.matrix(norm_sub))))
    scores   <- colMeans(z_mat, na.rm=TRUE)

    write.csv(data.frame(score=scores, row.names=names(scores)),
              file.path(out$csv, "signature_values", paste0(sig_name, ".csv")))

    col_name <- paste0("sig_", sig_name)
    colData(cds)[[col_name]] <- scores[colnames(cds)]

    ggsave(file.path(out$plots, paste0("umap_", sig_name, ".pdf")),
      plot_cells(cds, color_cells_by=col_name, show_trajectory_graph=FALSE,
                 label_cell_groups=FALSE) +
        scale_color_viridis_c(option="plasma", name="Score") +
        ggtitle(paste("Signature:", sig_name)),
      width=8, height=6)

    if ("gene_comb" %in% names(colData(cds))) {
      summary_rows[[sig_name]] <- data.frame(
        gene_comb=colData(cds)$gene_comb, score=scores[colnames(cds)]) %>%
        dplyr::group_by(gene_comb) %>%
        dplyr::summarise(mean_score=mean(score,na.rm=TRUE), n_cells=n(), .groups="drop") %>%
        dplyr::mutate(signature=sig_name)
    }
  }

  writeLines(coverage_log, file.path(out$stats, "signature_coverage.txt"))

  if (length(summary_rows) > 0) {
    df <- dplyr::bind_rows(summary_rows)
    write.csv(df, file.path(out$csv, "signature_summary.csv"), row.names=FALSE)
    log_info("Signature summary saved.")
  }

  log_info("Step 08 complete.")
  log_info("Output folder: ", out$out_dir)
}

main()
