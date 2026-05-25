#!/usr/bin/env Rscript
# =============================================================================
# CIRI Pipeline — Step 06: Subclustering & Trajectory
# =============================================================================
# INPUT  (auto-resolved from Step 03 to_scratch/):
#   processed_cds.RData
#
# OUTPUT STRUCTURE:
#   Output/
#     <YYMMDD>_step06_trajectory_<sample>/
#       plots/      UMAP_<group>_mainclusters.pdf
#                   UMAP_<group>_subclusters.pdf
#                   UMAP_<group>_trajectory.pdf
#                   UMAP_<group>_pseudotime.pdf
#       csv/        pseudotime_<group>.csv
#       R_objects/  processed_cds_<group>.RData
#       to_scratch/ processed_cds_<group>.RData   ← consumed by Steps 07–08
#                   pseudotime_<group>.csv         ← consumed by Step 07
#
# Usage:
#   Rscript steps/step_06_subclusters_trajectory.R \
#     --sample     AB011 \
#     --clusters   5 \
#     --root_gene  SOX2 \
#     --group      muscle
#   Rscript steps/step_06_subclusters_trajectory.R --help
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
                     help="Cluster(s) to subset, dash-separated: '5' or '3-4'"),
  root_gene   = list(flag="--root_gene",   required=TRUE,  type="character",
                     help="Gene marking the trajectory root (highest-expression node)"),
  group       = list(flag="--group",       required=TRUE,  type="character",
                     help="Short name for this lineage, used in filenames"),
  resolution  = list(flag="--resolution",  required=FALSE, type="numeric",
                     default=1e-3,      help="Re-clustering resolution within the subset"),
  n_dims      = list(flag="--n_dims",      required=FALSE, type="integer",
                     default=50L,       help="PCA dimensions for subset"),
  seed        = list(flag="--seed",        required=FALSE, type="integer",
                     default=42L,       help="Random seed")
)

main <- function() {
  suppressMessages({
    library(monocle3); library(dplyr); library(ggplot2); library(viridis)
  })

  p   <- parse_args(SPEC)
  out <- make_out_dirs(p$output_root, "step06_trajectory", p$sample)
  set.seed(p$seed)

  scratch03       <- if (!is.null(p$scratch)) p$scratch else
    find_scratch(p$output_root, "step03_load", p$sample)
  cds_path        <- file.path(scratch03, "processed_cds.RData")
  clusters_to_keep <- unlist(strsplit(p$clusters, "-"))

  step_banner("06", "Subclustering & Trajectory",
    inputs  = c(cds_path),
    outputs = c(out$plots, out$csv, out$R_objects,
                file.path(out$to_scratch, paste0("processed_cds_", p$group, ".RData")),
                file.path(out$to_scratch, paste0("pseudotime_", p$group, ".csv")))
  )
  assert_file(cds_path, hint = "Run step_03_load.R first.")

  cds <- load_rdata(cds_path)
  colData(cds)$clusters <- clusters(cds)
  log_info("Subsetting cluster(s): ", paste(clusters_to_keep, collapse=", "))

  cds_sub <- cds[, colData(cds)$clusters %in% clusters_to_keep]
  colData(cds_sub)$clusters_main <- colData(cds_sub)$clusters
  log_info("Cells in subset: ", ncol(cds_sub))

  log_info("Re-preprocessing subset ...")
  cds_sub <- preprocess_cds(cds_sub, num_dim=p$n_dims)
  cds_sub <- reduce_dimension(cds_sub, reduction_method="UMAP")
  cds_sub <- cluster_cells(cds_sub, resolution=p$resolution, random_seed=p$seed)
  colData(cds_sub)$clusters_sub <- clusters(cds_sub)
  log_info("Sub-clusters: ", length(unique(clusters(cds_sub))))

  save_umap_sub <- function(color_by, fname, label=FALSE, traj=FALSE) {
    pp <- plot_cells(cds_sub, color_cells_by=color_by, show_trajectory_graph=traj,
                     label_cell_groups=label, label_leaves=label, label_branch_points=label) +
      ggtitle(paste("UMAP —", color_by, "|", p$group))
    ggsave(file.path(out$plots, fname), pp, width=6, height=5)
  }
  save_umap_sub("clusters_main", paste0("UMAP_", p$group, "_mainclusters.pdf"),  label=TRUE)
  save_umap_sub("clusters_sub",  paste0("UMAP_", p$group, "_subclusters.pdf"),   label=TRUE)

  log_info("Learning trajectory graph ...")
  cds_sub <- learn_graph(cds_sub)

  # Set root by root_gene expression
  root_node <- NULL
  if (p$root_gene %in% rownames(cds_sub)) {
    root_expr <- as.numeric(normalized_counts(cds_sub)[p$root_gene, ])
    graph_nodes <- igraph::V(principal_graph(cds_sub)[["UMAP"]])$name
    closest     <- cds_sub@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
    mean_per_node <- sapply(graph_nodes, function(nd) {
      idx <- which(closest == nd)
      if (!length(idx)) return(-Inf)
      mean(root_expr[idx])
    })
    root_node <- graph_nodes[which.max(mean_per_node)]
    log_info("Root node (", p$root_gene, "): ", root_node)
  } else log_warn("root_gene '", p$root_gene, "' not in matrix; root will be arbitrary.")

  cds_sub <- tryCatch(
    order_cells(cds_sub, root_pr_nodes=root_node),
    error = function(e) { log_warn("order_cells failed: ", conditionMessage(e), " — running without root."); order_cells(cds_sub) }
  )

  save_umap_sub("pseudotime",   paste0("UMAP_", p$group, "_pseudotime.pdf"),  traj=TRUE)
  save_umap_sub("clusters_sub", paste0("UMAP_", p$group, "_trajectory.pdf"),  traj=TRUE, label=TRUE)

  # csv/ + to_scratch/
  pt_df   <- data.frame(pseudotime=pseudotime(cds_sub), row.names=colnames(cds_sub))
  pt_file <- paste0("pseudotime_", p$group, ".csv")
  write.csv(pt_df, file.path(out$csv, pt_file))
  file.copy(file.path(out$csv, pt_file), file.path(out$to_scratch, pt_file))

  # R_objects/ + to_scratch/
  cds_file <- paste0("processed_cds_", p$group, ".RData")
  rdata_path <- file.path(out$R_objects, cds_file)
  save(cds_sub, file=rdata_path)
  file.copy(rdata_path, file.path(out$to_scratch, cds_file))

  log_info("Step 06 complete.")
  log_info("Output folder: ", out$out_dir)
}

main()
