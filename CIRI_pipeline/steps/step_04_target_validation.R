#!/usr/bin/env Rscript
# =============================================================================
# CIRI Pipeline — Step 04: Target Validation
# =============================================================================
# Checks knockdown / activation efficiency per targeted gene vs control.
#
# INPUT  (auto-resolved from Step 03 to_scratch/):
#   processed_cds.RData
#
# OUTPUT STRUCTURE:
#   Output/
#     <YYMMDD>_step04_validation_<sample>/
#       plots/   validation_violin_<gene>.pdf
#       csv/     validation_summary.csv
#       stats/   validation_stats.txt
#
# Usage:
#   Rscript steps/step_04_target_validation.R \
#     --sample    AB011 \
#     --control_a "NTCa-NA"
#   Rscript steps/step_04_target_validation.R --help
# =============================================================================

STEPS_DIR <- dirname(sys.frame(1)$ofile)
source(file.path(STEPS_DIR, "utils.R"))

SPEC <- list(
  output_root = list(flag="--output_root", required=FALSE, type="character",
                     default="Output",  help="Top-level output directory"),
  sample      = list(flag="--sample",      required=FALSE, type="character",
                     default="CIRI",    help="Sample / experiment name"),
  scratch     = list(flag="--scratch",     required=FALSE, type="character",
                     default=NULL,      help="Path to Step 03 to_scratch/ (auto-resolved if omitted)"),
  control_a   = list(flag="--control_a",   required=FALSE, type="character",
                     default="NTCa-NA", help="Control gene_comb string for CRISPRa"),
  control_i   = list(flag="--control_i",   required=FALSE, type="character",
                     default="NTCa-NA", help="Control gene_comb string for CRISPRi"),
  min_cells   = list(flag="--min_cells",   required=FALSE, type="integer",
                     default=5L,        help="Min cells per group to plot")
)

main <- function() {
  suppressMessages({
    library(monocle3); library(dplyr); library(tidyr)
    library(ggplot2); library(ggsignif); library(viridis)
  })

  p   <- parse_args(SPEC)
  out <- make_out_dirs(p$output_root, "step04_validation", p$sample)

  scratch03 <- if (!is.null(p$scratch)) p$scratch else
    find_scratch(p$output_root, "step03_load", p$sample)
  cds_path  <- file.path(scratch03, "processed_cds.RData")

  step_banner("04", "Target Validation",
    inputs  = c(cds_path),
    outputs = c(out$plots, out$csv, out$stats)
  )
  assert_file(cds_path, hint = "Run step_03_load.R first.")

  cds      <- load_rdata(cds_path)
  norm_mat <- normalized_counts(cds)

  cell_df <- data.frame(nomi=colnames(norm_mat), row.names=colnames(norm_mat)) %>%
    tidyr::separate(nomi, into=c("cellID","sample","guide_a","guide_i"),
                    sep="-", remove=FALSE, fill="right") %>%
    dplyr::mutate(
      guide_a   = na_if(guide_a, "NA"),
      guide_i   = na_if(guide_i, "NA"),
      gene_a    = sapply(strsplit(guide_a,"_"), `[`, 1),
      gene_i    = sapply(strsplit(guide_i,"_"), `[`, 1),
      gene_comb = paste(gene_a, gene_i, sep="-")
    )

  target_a <- unique(cell_df$gene_a[!is.na(cell_df$gene_a) & cell_df$gene_comb != p$control_a])
  target_i <- unique(cell_df$gene_i[!is.na(cell_df$gene_i) & cell_df$gene_comb != p$control_i])
  gene_class <- dplyr::bind_rows(
    data.frame(gene=target_a, class="CRISPRa"),
    data.frame(gene=target_i, class="CRISPRi")
  ) %>% dplyr::filter(!is.na(gene)) %>% dplyr::distinct()

  log_info("CRISPRa targets: ", paste(target_a, collapse=", "))
  log_info("CRISPRi targets: ", paste(target_i, collapse=", "))

  summary_rows <- list()
  for (i in seq_len(nrow(gene_class))) {
    gene  <- gene_class$gene[i]
    class <- gene_class$class[i]
    if (!gene %in% rownames(norm_mat)) { log_warn("Gene not in matrix: ", gene); next }
    expr_vec <- norm_mat[gene, ]
    perturbed <- if (class=="CRISPRa") cell_df$nomi[!is.na(cell_df$gene_a) & cell_df$gene_a==gene] \
                 else                  cell_df$nomi[!is.na(cell_df$gene_i) & cell_df$gene_i==gene]
    control   <- if (class=="CRISPRa") cell_df$nomi[cell_df$gene_comb==p$control_a] \
                 else                  cell_df$nomi[cell_df$gene_comb==p$control_i]
    if (length(perturbed) < p$min_cells || length(control) < p$min_cells) {
      log_warn("Too few cells for ", gene, " — skipped"); next
    }
    plot_df <- dplyr::bind_rows(
      data.frame(expr=as.numeric(expr_vec[perturbed]), group=paste0(gene," perturbed")),
      data.frame(expr=as.numeric(expr_vec[control]),   group="Control")
    )
    ggsave(file.path(out$plots, paste0("validation_violin_", gene, ".pdf")),
      ggplot(plot_df, aes(x=group, y=expr, fill=group)) +
        geom_violin(trim=FALSE, alpha=0.7) +
        geom_boxplot(width=0.1, fill="white", alpha=0.7, outlier.size=0.5) +
        theme_minimal(base_size=13) + labs(title=paste(class,"—",gene), x=NULL, y="Norm. expression") +
        theme(legend.position="none"),
      width=5, height=5)
    summary_rows[[gene]] <- data.frame(
      gene=gene, class=class, n_perturbed=length(perturbed), n_control=length(control),
      median_perturbed=median(as.numeric(expr_vec[perturbed])),
      median_control  =median(as.numeric(expr_vec[control])),
      fold_change     =median(as.numeric(expr_vec[perturbed])) /
                       (median(as.numeric(expr_vec[control]))+1e-9))
  }

  if (length(summary_rows) > 0) {
    df <- dplyr::bind_rows(summary_rows)
    write.csv(df, file.path(out$csv, "validation_summary.csv"), row.names=FALSE)
    writeLines(capture.output(print(df)), file.path(out$stats, "validation_stats.txt"))
    log_info("Validation summary:"); print(df)
  }

  log_info("Step 04 complete.")
  log_info("Output folder: ", out$out_dir)
}

main()
