#!/usr/bin/env Rscript
# =============================================================================
# CIRI Pipeline — Step 03: Load & Preprocess
# =============================================================================
# INPUT  (auto-resolved from Step 02 to_scratch/):
#   annotated_matrix.csv
#
# OUTPUT STRUCTURE:
#   Output/
#     <YYMMDD>_step03_load_<sample>/
#       plots/      umap_clusters.pdf, umap_sample.pdf, umap_guide_a/i.pdf
#       stats/      cluster_summary.txt
#       R_objects/  processed_cds.RData
#       to_scratch/ processed_cds.RData   ← consumed by Steps 04–06
#
# Usage:
#   Rscript steps/step_03_load.R \
#     --sample      AB011 \
#     --resolution  5e-5
#   Rscript steps/step_03_load.R --help
# =============================================================================

STEPS_DIR <- dirname(sys.frame(1)$ofile)
source(file.path(STEPS_DIR, "utils.R"))

SPEC <- list(
  output_root = list(flag="--output_root", required=FALSE, type="character",
                     default="Output",   help="Top-level output directory"),
  sample      = list(flag="--sample",      required=FALSE, type="character",
                     default="CIRI",     help="Sample / experiment name"),
  scratch     = list(flag="--scratch",     required=FALSE, type="character",
                     default=NULL,
                     help="Path to Step 02 to_scratch/ (auto-resolved if omitted)"),
  resolution  = list(flag="--resolution",  required=FALSE, type="numeric",
                     default=5e-5,
                     help="Monocle3 Leiden clustering resolution (try 1e-5 → 1e-3)"),
  n_dims      = list(flag="--n_dims",      required=FALSE, type="integer",
                     default=100L,       help="PCA dimensions"),
  seed        = list(flag="--seed",        required=FALSE, type="integer",
                     default=1234597698L, help="Random seed")
)

main <- function() {
  suppressMessages({
    library(monocle3); library(dplyr); library(tidyr)
    library(ggplot2); library(viridis); library(ggrepel); library(gtools)
  })

  p   <- parse_args(SPEC)
  out <- make_out_dirs(p$output_root, "step03_load", p$sample)
  set.seed(p$seed)

  scratch02 <- if (!is.null(p$scratch)) p$scratch else
    find_scratch(p$output_root, "step02_filter", p$sample)
  mat_path  <- file.path(scratch02, "annotated_matrix.csv")

  step_banner("03", "Load & Preprocess",
    inputs  = c(mat_path),
    outputs = c(out$plots, out$R_objects,
                file.path(out$to_scratch, "processed_cds.RData"))
  )
  assert_file(mat_path, hint = "Run step_02_anno_filter.R first.")

  log_info("Reading annotated matrix ...")
  exp           <- read.csv(mat_path, header=TRUE, check.names=FALSE)
  rownames(exp) <- exp[,1]; exp[,1] <- NULL
  log_info("Matrix: ", nrow(exp), " genes × ", ncol(exp), " cells")

  cell_meta <- data.frame(nomi=colnames(exp), row.names=colnames(exp)) %>%
    tidyr::separate(nomi, into=c("cellID","sample","guide_a","guide_i"),
                    sep="-", remove=FALSE, fill="right") %>%
    dplyr::mutate(
      guide_a   = na_if(guide_a, "NA"),
      guide_i   = na_if(guide_i, "NA"),
      comb      = paste(guide_a, guide_i, sep="-"),
      gene_a    = sapply(strsplit(guide_a,"_"), `[`, 1),
      gene_i    = sapply(strsplit(guide_i,"_"), `[`, 1),
      gene_comb = paste(gene_a, gene_i, sep="-")
    )

  gene_meta <- data.frame(gene_short_name=rownames(exp), row.names=rownames(exp))

  log_info("Creating Monocle3 cell_data_set ...")
  cds <- new_cell_data_set(as.matrix(exp), cell_metadata=cell_meta, gene_metadata=gene_meta)

  log_info("Preprocessing (PCA, ", p$n_dims, " dims) ...")
  cds <- preprocess_cds(cds, num_dim=p$n_dims)

  log_info("Running UMAP ...")
  cds <- reduce_dimension(cds)

  log_info("Clustering (resolution=", p$resolution, ") ...")
  cds <- cluster_cells(cds, resolution=p$resolution)
  colData(cds)$clusters  <- clusters(cds)
  n_clust <- length(unique(clusters(cds)))
  log_info("Found ", n_clust, " clusters")

  # plots/
  save_umap <- function(color_by, fname, label=FALSE) {
    pp <- plot_cells(cds, color_cells_by=color_by, show_trajectory_graph=FALSE,
                     label_cell_groups=label) + ggtitle(paste("UMAP —", color_by))
    ggsave(file.path(out$plots, fname), pp, width=10, height=8)
  }
  save_umap("cluster",  "umap_clusters.pdf", label=TRUE)
  save_umap("sample",   "umap_sample.pdf")
  if ("gene_a" %in% names(colData(cds))) save_umap("gene_a", "umap_guide_a.pdf")
  if ("gene_i" %in% names(colData(cds))) save_umap("gene_i", "umap_guide_i.pdf")

  # stats/
  clust_tbl <- as.data.frame(table(clusters(cds)))
  names(clust_tbl) <- c("cluster","n_cells")
  write.csv(clust_tbl, file.path(out$stats, "cluster_sizes.csv"), row.names=FALSE)
  writeLines(c(paste("Resolution:", p$resolution), paste("N clusters:", n_clust),
               capture.output(print(clust_tbl))),
             file.path(out$stats, "cluster_summary.txt"))

  # R_objects/ + to_scratch/
  rdata_path <- file.path(out$R_objects, "processed_cds.RData")
  save(cds, file=rdata_path)
  file.copy(rdata_path, file.path(out$to_scratch, "processed_cds.RData"))

  log_info("Step 03 complete.")
  log_info("Output folder: ", out$out_dir)
}

main()
