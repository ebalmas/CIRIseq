#!/usr/bin/env Rscript
# =============================================================================
# CIRI Pipeline — Step 07: Pseudotime Statistics (KS Test)
# =============================================================================
# INPUT  (auto-resolved from Step 06 to_scratch/):
#   processed_cds_<group>.RData
#   pseudotime_<group>.csv
#
# OUTPUT STRUCTURE:
#   Output/
#     <YYMMDD>_step07_pseudotime_<sample>/
#       csv/    ks_results_<group>.csv
#               ks_results_<group>_<sample>.csv  (if --run_per_sample)
#       plots/  volcano_<group>.pdf
#               ecdf_<hit>_<group>.pdf
#       stats/  ks_summary_<group>.txt
#
# Usage:
#   Rscript steps/step_07_pseudotime_stats.R \
#     --sample   AB011 \
#     --group    muscle \
#     --control  "NTCa-NA"
#   Rscript steps/step_07_pseudotime_stats.R --help
# =============================================================================

STEPS_DIR <- dirname(sys.frame(1)$ofile)
source(file.path(STEPS_DIR, "utils.R"))

SPEC <- list(
  output_root    = list(flag="--output_root",    required=FALSE, type="character",
                        default="Output",   help="Top-level output directory"),
  sample         = list(flag="--sample",         required=FALSE, type="character",
                        default="CIRI",     help="Sample / experiment name"),
  scratch        = list(flag="--scratch",        required=FALSE, type="character",
                        default=NULL,       help="Path to Step 06 to_scratch/ (auto-resolved if omitted)"),
  group          = list(flag="--group",          required=TRUE,  type="character",
                        help="Group name matching Step 06 --group"),
  control        = list(flag="--control",        required=TRUE,  type="character",
                        help="Control gene_comb string, e.g. 'NTCa-NA'"),
  min_cells      = list(flag="--min_cells",      required=FALSE, type="integer",
                        default=8L,         help="Min cells per perturbation for KS test"),
  run_per_sample = list(flag="--run_per_sample", required=FALSE, type="logical",
                        default=FALSE,      help="Also run per-sample analysis"),
  ecdf_top_n     = list(flag="--ecdf_top_n",    required=FALSE, type="integer",
                        default=10L,        help="Top N hits to get individual ECDF plots")
)

run_ks <- function(pt_df, control_label, min_cells) {
  ctrl_pt <- pt_df$pseudotime[pt_df$gene_comb == control_label]
  ctrl_pt <- ctrl_pt[!is.na(ctrl_pt) & is.finite(ctrl_pt)]
  if (length(ctrl_pt) < min_cells)
    stop(paste("Control '", control_label, "' has < ", min_cells, " finite pseudotime values."))
  rows <- lapply(unique(pt_df$gene_comb), function(comb) {
    if (comb == control_label) return(NULL)
    pert_pt <- pt_df$pseudotime[pt_df$gene_comb == comb]
    pert_pt <- pert_pt[!is.na(pert_pt) & is.finite(pert_pt)]
    if (length(pert_pt) < min_cells) return(NULL)
    ks <- ks.test(pert_pt, ctrl_pt)
    data.frame(gene_comb=comb, n_cells=length(pert_pt), ks_stat=ks$statistic,
               p_value=ks$p.value,
               direction=ifelse(median(pert_pt)>median(ctrl_pt),"faster","slower"),
               median_pert=median(pert_pt), median_ctrl=median(ctrl_pt))
  })
  df <- dplyr::bind_rows(rows)
  df$p_adj <- p.adjust(df$p_value, method="BH")
  dplyr::arrange(df, p_adj, desc(ks_stat))
}

main <- function() {
  suppressMessages({
    library(monocle3); library(dplyr); library(ggplot2); library(ggrepel); library(viridis)
  })

  p   <- parse_args(SPEC)
  out <- make_out_dirs(p$output_root, "step07_pseudotime", p$sample)

  scratch06 <- if (!is.null(p$scratch)) p$scratch else
    find_scratch(p$output_root, "step06_trajectory", p$sample)
  cds_path <- file.path(scratch06, paste0("processed_cds_", p$group, ".RData"))
  pt_path  <- file.path(scratch06, paste0("pseudotime_", p$group, ".csv"))

  step_banner("07", "Pseudotime Statistics (KS Test)",
    inputs  = c(cds_path, pt_path),
    outputs = c(out$csv, out$plots, out$stats)
  )
  assert_file(cds_path, hint = "Run step_06_subclusters_trajectory.R first.")
  assert_file(pt_path,  hint = "Run step_06_subclusters_trajectory.R first.")

  cds     <- load_rdata(cds_path)
  pt      <- read.csv(pt_path, row.names=1); names(pt) <- "pseudotime"
  cell_df <- as.data.frame(colData(cds))
  shared  <- intersect(rownames(pt), rownames(cell_df))
  if (!length(shared)) stop("No matching cell names between CDS and pseudotime CSV.")
  cell_df <- cell_df[shared, ]
  cell_df$pseudotime <- pt[shared, "pseudotime"]

  log_info("Running KS tests (control='", p$control, "', min_cells=", p$min_cells, ") ...")
  ks_res  <- run_ks(cell_df, p$control, p$min_cells)
  out_csv <- file.path(out$csv, paste0("ks_results_", p$group, ".csv"))
  write.csv(ks_res, out_csv, row.names=FALSE)
  n_sig   <- sum(ks_res$p_adj < 0.05, na.rm=TRUE)
  log_info("Significant hits (p.adj < 0.05): ", n_sig)

  # Volcano
  ks_plot <- dplyr::mutate(ks_res,
    log10_padj  = -log10(p_adj + 1e-300),
    label       = ifelse(p_adj < 0.05 & dplyr::row_number() <= 20, gene_comb, NA_character_))
  ggsave(file.path(out$plots, paste0("volcano_", p$group, ".pdf")),
    ggplot(ks_plot, aes(x=ks_stat, y=log10_padj, colour=direction, label=label)) +
      geom_point(aes(size=n_cells), alpha=0.7) +
      ggrepel::geom_text_repel(na.rm=TRUE, size=3, max.overlaps=20) +
      geom_hline(yintercept=-log10(0.05), linetype="dashed", colour="grey50") +
      scale_colour_manual(values=c(faster="#2166AC", slower="#D6604D")) +
      scale_size_continuous(range=c(1,5)) + theme_minimal(base_size=13) +
      labs(title=paste("Pseudotime KS test —", p$group), x="KS statistic",
           y="-log10(adj. p-value)", colour="Direction", size="n cells"),
    width=10, height=8)

  # ECDF top N
  ctrl_pt <- cell_df$pseudotime[cell_df$gene_comb==p$control]
  ctrl_pt <- ctrl_pt[!is.na(ctrl_pt) & is.finite(ctrl_pt)]
  for (comb in head(ks_res$gene_comb, p$ecdf_top_n)) {
    pert_pt <- cell_df$pseudotime[cell_df$gene_comb==comb]
    pert_pt <- pert_pt[!is.na(pert_pt) & is.finite(pert_pt)]
    ecdf_df <- dplyr::bind_rows(
      data.frame(pseudotime=pert_pt, group=comb),
      data.frame(pseudotime=ctrl_pt, group=p$control))
    safe_name <- gsub("[^A-Za-z0-9_.-]","_", comb)
    ggsave(file.path(out$plots, paste0("ecdf_", safe_name, "_", p$group, ".pdf")),
      ggplot(ecdf_df, aes(x=pseudotime, colour=group)) +
        stat_ecdf(linewidth=1) + theme_minimal(base_size=12) +
        scale_colour_manual(values=c("steelblue","grey50"), breaks=c(comb,p$control)) +
        labs(title=paste("ECDF —", comb, "vs", p$control), x="Pseudotime",
             y="Cumulative fraction", colour=NULL),
      width=7, height=5)
  }

  if (p$run_per_sample && "sample" %in% names(cell_df)) {
    for (samp in unique(cell_df$sample)) {
      tryCatch({
        ks_s <- run_ks(dplyr::filter(cell_df, sample==samp), p$control, p$min_cells)
        write.csv(ks_s, file.path(out$csv, paste0("ks_results_", p$group, "_", samp, ".csv")),
                  row.names=FALSE)
        log_info("Per-sample KS saved: ", samp)
      }, error = function(e) log_warn("  Sample ", samp, " skipped: ", conditionMessage(e)))
    }
  }

  writeLines(c(paste("Group:", p$group), paste("Control:", p$control),
               paste("N tested:", nrow(ks_res)), paste("N significant (p.adj<0.05):", n_sig),
               "Top 10:", capture.output(print(head(ks_res,10)))),
             file.path(out$stats, paste0("ks_summary_", p$group, ".txt")))

  log_info("Step 07 complete.")
  log_info("Output folder: ", out$out_dir)
}

main()
