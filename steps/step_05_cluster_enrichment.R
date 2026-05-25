#!/usr/bin/env Rscript
# =============================================================================
# CIRI Pipeline — Step 05: Cluster Enrichment Analysis
# =============================================================================
# INPUT  (auto-resolved from Step 03 to_scratch/):
#   processed_cds.RData
#
# OUTPUT STRUCTURE:
#   Output/
#     <YYMMDD>_step05_enrichment_<sample>/
#       csv/    cluster_enrichment_<group>.csv
#       plots/  gene_a_in_clusters.pdf, gene_i_in_clusters.pdf
#               heatmap_enrichment_<group>.pdf
#               scatter_enrichment_<group>.pdf
#       stats/  enrichment_summary.txt
#
# Usage:
#   Rscript steps/step_05_cluster_enrichment.R \
#     --sample   AB011 \
#     --clusters 5 \
#     --control  "NTCa-NA"
#   Rscript steps/step_05_cluster_enrichment.R --help
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
  clusters    = list(flag="--clusters",    required=TRUE,  type="character",
                     help="Target cluster ID(s), dash-separated: '5' or '3-4'"),
  control     = list(flag="--control",     required=TRUE,  type="character",
                     help="Control gene_comb label, e.g. 'NTCa-NA'"),
  min_cells   = list(flag="--min_cells",   required=FALSE, type="integer",
                     default=10L,       help="Min cells per perturbation to include")
)

main <- function() {
  suppressMessages({
    library(monocle3); library(dplyr); library(tidyr)
    library(ggplot2); library(viridis); library(ggrepel)
  })

  p   <- parse_args(SPEC)
  out <- make_out_dirs(p$output_root, "step05_enrichment", p$sample)

  scratch03       <- if (!is.null(p$scratch)) p$scratch else
    find_scratch(p$output_root, "step03_load", p$sample)
  cds_path        <- file.path(scratch03, "processed_cds.RData")
  target_clusters <- unlist(strsplit(p$clusters, "-"))
  group_label     <- paste("Group", paste(target_clusters, collapse="_"), sep="_")

  step_banner("05", "Cluster Enrichment Analysis",
    inputs  = c(cds_path),
    outputs = c(out$csv, out$plots, out$stats)
  )
  assert_file(cds_path, hint = "Run step_03_load.R first.")

  cds <- load_rdata(cds_path)
  colData(cds)$clusters <- clusters(cds)
  df  <- as.data.frame(colData(cds))

  # Overview bar charts
  overview_plot <- function(df, gene_col, fname) {
    if (!gene_col %in% names(df)) return(invisible(NULL))
    counts <- df %>% group_by(sample, clusters, .data[[gene_col]]) %>%
      summarise(count=n(), .groups="drop") %>%
      group_by(sample, clusters) %>%
      mutate(total=sum(count), pct=count/total*100) %>% ungroup()
    n_vals <- length(unique(counts[[gene_col]]))
    p <- ggplot(counts, aes(x=factor(sample), y=pct, fill=.data[[gene_col]])) +
      geom_bar(stat="identity", colour="white", linewidth=0.2) +
      scale_fill_manual(values=viridis::turbo(n_vals)) +
      facet_grid(~clusters) + theme_minimal() +
      labs(title=paste(gene_col,"across clusters"), x="Sample", y="% cells") +
      theme(axis.text.x=element_text(angle=45, hjust=1))
    ggsave(file.path(out$plots, fname), p, width=20, height=10)
  }
  overview_plot(df, "gene_a", "gene_a_in_clusters.pdf")
  overview_plot(df, "gene_i", "gene_i_in_clusters.pdf")

  # Enrichment in target clusters
  enrichment <- df %>%
    mutate(in_target = clusters %in% target_clusters) %>%
    group_by(sample, gene_comb) %>%
    summarise(n_total=n(), n_target=sum(in_target),
              pct_target=n_target/n_total*100, .groups="drop") %>%
    filter(n_total >= p$min_cells)

  log_info("Perturbations with >= ", p$min_cells, " cells: ", length(unique(enrichment$gene_comb)))

  write.csv(enrichment, file.path(out$csv, paste0("cluster_enrichment_", group_label, ".csv")),
            row.names=FALSE)

  # Heatmap
  samples <- unique(enrichment$sample)
  if (length(samples) >= 2) {
    wide <- enrichment %>% select(sample, gene_comb, pct_target) %>%
      tidyr::pivot_wider(names_from=sample, values_from=pct_target, values_fill=NA)
    complete_comb <- wide %>%
      filter(rowSums(is.na(across(-gene_comb)))==0) %>% pull(gene_comb)
    heat_df <- filter(enrichment, gene_comb %in% complete_comb)
    ggsave(file.path(out$plots, paste0("heatmap_enrichment_", group_label, ".pdf")),
      ggplot(heat_df, aes(x=factor(sample), y=gene_comb, fill=pct_target)) +
        geom_tile(colour="white", linewidth=0.3) +
        scale_fill_viridis_c(option="plasma", name="% in target") +
        theme_minimal(base_size=10) +
        labs(title=paste("Enrichment — cluster(s)", paste(target_clusters, collapse="+")),
             x="Sample", y="Perturbation"),
      width=8, height=max(5, length(complete_comb)*0.25))

    if (length(samples)==2) {
      scatter_df <- wide %>% filter(gene_comb %in% complete_comb)
      names(scatter_df)[2:3] <- c("s1","s2")
      ctrl_row <- filter(scatter_df, gene_comb==p$control)
      label_df <- filter(scatter_df, gene_comb!=p$control) %>%
        arrange(desc(abs(s1-s2))) %>% head(20)
      ggsave(file.path(out$plots, paste0("scatter_enrichment_", group_label, ".pdf")),
        ggplot(scatter_df, aes(x=s1, y=s2, label=gene_comb)) +
          geom_point(alpha=0.5, colour="steelblue") +
          geom_point(data=ctrl_row, colour="red", size=3) +
          ggrepel::geom_text_repel(data=label_df, size=3, max.overlaps=20) +
          geom_abline(linetype="dashed", colour="grey60") + theme_minimal() +
          labs(title=paste("Replicate comparison — cluster(s)", paste(target_clusters, collapse="+")),
               x=paste("Sample", samples[1]), y=paste("Sample", samples[2])),
        width=8, height=8)
    }
  }

  ctrl_pct <- mean(filter(enrichment, gene_comb==p$control)$pct_target, na.rm=TRUE)
  top5     <- enrichment %>% group_by(gene_comb) %>%
    summarise(mean_pct=mean(pct_target), .groups="drop") %>%
    arrange(desc(mean_pct)) %>% head(5)
  writeLines(c(paste("Control mean % in target:", round(ctrl_pct,2)),
               "Top 5 perturbations:", capture.output(print(top5))),
             file.path(out$stats, "enrichment_summary.txt"))

  log_info("Step 05 complete.")
  log_info("Output folder: ", out$out_dir)
}

main()
