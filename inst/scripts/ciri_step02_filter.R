#!/usr/bin/env Rscript
# =============================================================================
# ciri_step02_filter.R
# =============================================================================
# Filters the Seurat object from ciri_step02_annotate.R, then builds a
# Monocle3 CDS at multiple clustering resolutions.
#
# RUN THIS SCRIPT THREE TIMES:
# ─────────────────────────────────────────────────────────────────────────────
# RUN 1 — filter + build Monocle3 + see all 5 resolutions
#
#   Rscript ciri_step02_filter.R \
#     --sample    AB014_AB016 \
#     --mito_hi   10  --mito_lo 0  --ribo_lo 1 \
#     --nGene_lo  300 --nGene_hi 7000 --nUMI_lo 100
#
#   Then inspect:
#     → Output/monocle/<date>/umap/<date>_variance_knee_plot_dim30.png
#       Does the curve flatten before dim 30? If yes, 30 is fine.
#       If it keeps dropping, re-run with --num_dim 50 (or higher).
#
# RUN 2 (only if knee plot says num_dim is wrong)
#
#   Rscript ciri_step02_filter.R --sample AB014_AB016 --num_dim 50 [+ same thresholds]
#
#   Then inspect:
#     → Output/monocle/<date>/umap/<date>_cds_1 through cds_5 PDFs
#       Compare the 5 UMAP resolutions and pick the one that best
#       separates your biology (e.g. cds_3).
#
# RUN 3 — finalise the chosen resolution
#
#   Rscript ciri_step02_filter.R --sample AB014_AB016 --chosen_cds cds_3 [+ same args]
#
#   This saves cds_final.rds to to_scratch/ for downstream steps.
# ─────────────────────────────────────────────────────────────────────────────
#
# ALL PARAMETERS can be passed as command-line arguments.
# You never need to edit this script — just change the arguments.
#
# FULL ARGUMENT LIST:
#   --sample       experiment name (required)
#   --scratch_dir  folder with seurat_annotated.rds (default: scratch)
#   --output_root  top-level output folder (default: Output)
#   --mito_lo      percent.mt >= X      (default: 0)
#   --mito_hi      percent.mt <  X      (default: 10)
#   --ribo_lo      percent.ribo >= X    (default: 1)
#   --nGene_lo     nFeature_RNA > X     (default: 300)
#   --nGene_hi     nFeature_RNA < X     (default: 7000)
#   --nUMI_lo      nCount_RNA > X       (default: 100)
#   --num_dim      PCA dimensions       (default: 30)
#   --chosen_cds   cds_1 … cds_7        (default: cds_3)
#                  only matters on Run 3 — all 7 are always computed and saved
#
# OUTPUT STRUCTURE
# ─────────────────────────────────────────────────────────────────────────────
#   Output/QC/<date>/
#     csv/        <date>_post_QCcounts.csv
#                 <date>_cell_metadata_filtered.csv
#     filtering/  <date>_all_QC_after_FILTER.pdf
#     R_objects/  <date>_seurat_data_filtered.RData
#   Output/monocle/<date>/
#     umap/       <date>_variance_knee_plot_dim<N>.png    ← check after Run 1
#                 <date>_cds_1 … cds_5 UMAP PDFs         ← compare after Run 1/2
#                 <date>_UMAP.pdf                         ← chosen resolution
#     R_objects/  <date>_monocle_cds_clustered.RData
#                 <date>_monocle_cds_clustered_FINAL.RData
#     to_scratch/ cds_final.rds                          ← for downstream steps
# =============================================================================

# ---------------------------------------------------------------------------
# Defaults — all of these can be overridden from the command line.
# You should never need to edit this file directly.
# ---------------------------------------------------------------------------

# QC filter thresholds
FILTER_PARAMS <- list(
  mito_lo  = 0,
  mito_hi  = 10,
  ribo_lo  = 1,
  nGene_lo = 300,
  nGene_hi = 7000,
  nUMI_lo  = 100
)

# PCA dimensions for Monocle3
NUM_DIM <- 30

# Which clustering resolution to use as the final CDS
# (cds_1 = coarsest / fewest clusters … cds_5 = finest / most clusters)
# All 5 are always computed and saved — this only affects which one
# gets saved as cds_final.rds in to_scratch/
CHOSEN_CDS <- "cds_3"

# ---------------------------------------------------------------------------
# Misc config
# ---------------------------------------------------------------------------
INTERACTIVE_PARAMS <- list(
  scratch_dir  = "scratch",
  output_root  = "Output",
  sample       = "AB014_AB016"
)

pal <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
         "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
         "#920000","#924900","#db6d00","#24ff24","#ffff6d")

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
.parse_args <- function(defaults) {
  raw <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% raw || "-h" %in% raw) {
    cat("Rscript ciri_step02_filter.R\n",
        "  --scratch_dir <path>  where seurat_annotated.rds lives\n",
        "  --output_root <path>\n",
        "  --sample      <name>\n",
        "  --num_dim     <int>   PCA dims (default:", NUM_DIM, ")\n",
        "  --mito_lo/hi  --ribo_lo  --nGene_lo/hi  --nUMI_lo\n",
        "    override the FILTER_PARAMS at runtime\n")
    quit(status = 0)
  }
  p <- defaults; i <- 1
  while (i <= length(raw)) {
    switch(raw[i],
      "--scratch_dir" = { p$scratch_dir <- raw[i+1]; i <- i+2 },
      "--output_root" = { p$output_root <- raw[i+1]; i <- i+2 },
      "--sample"      = { p$sample      <- raw[i+1]; i <- i+2 },
      "--num_dim"     = { NUM_DIM       <<- as.integer(raw[i+1]); i <- i+2 },
      "--chosen_cds"  = { CHOSEN_CDS   <<- raw[i+1];             i <- i+2 },
      "--mito_lo"     = { FILTER_PARAMS$mito_lo  <<- as.numeric(raw[i+1]); i <- i+2 },
      "--mito_hi"     = { FILTER_PARAMS$mito_hi  <<- as.numeric(raw[i+1]); i <- i+2 },
      "--ribo_lo"     = { FILTER_PARAMS$ribo_lo  <<- as.numeric(raw[i+1]); i <- i+2 },
      "--nGene_lo"    = { FILTER_PARAMS$nGene_lo <<- as.numeric(raw[i+1]); i <- i+2 },
      "--nGene_hi"    = { FILTER_PARAMS$nGene_hi <<- as.numeric(raw[i+1]); i <- i+2 },
      "--nUMI_lo"     = { FILTER_PARAMS$nUMI_lo  <<- as.numeric(raw[i+1]); i <- i+2 },
      { i <- i+1 }
    )
  }
  p
}

p <- .parse_args(INTERACTIVE_PARAMS)
F <- FILTER_PARAMS  # short alias

# ---------------------------------------------------------------------------
# Libraries
# ---------------------------------------------------------------------------
suppressMessages({
  library(Seurat)
  library(monocle3)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(viridis)
})

theme_set(theme_bw(12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 15, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)))

log_info <- function(...) message("[INFO]  ", ...)

# ---------------------------------------------------------------------------
# Output dirs
# ---------------------------------------------------------------------------
d        <- format(Sys.Date(), "%y%m%d")
qc_root  <- file.path(p$output_root, "QC",      d)
mon_root <- file.path(p$output_root, "monocle",  d)
dirs <- list(
  qc_csv    = file.path(qc_root,  "csv"),
  qc_filter = file.path(qc_root,  "filtering"),
  qc_R      = file.path(qc_root,  "R_objects"),
  mon_R     = file.path(mon_root, "R_objects"),
  mon_umap  = file.path(mon_root, "umap"),
  mon_scr   = file.path(mon_root, "to_scratch")
)
for (dr in dirs) dir.create(dr, recursive = TRUE, showWarnings = FALSE)

pfx_qc_csv <- file.path(dirs$qc_csv,    d)
pfx_filter  <- file.path(dirs$qc_filter, d)
pfx_qc_R    <- file.path(dirs$qc_R,      d)
pfx_mon_R   <- file.path(dirs$mon_R,     d)
pfx_umap    <- file.path(dirs$mon_umap,  d)

sep <- strrep("=", 72)
message(sep)
message("  CIRI — Step 02b: Filter + Monocle3")
message(sprintf("  %s", Sys.time()))
message(sep)
message("  FILTER thresholds applied:")
message(sprintf("    percent.mt   : [%s, %s)", F$mito_lo,  F$mito_hi))
message(sprintf("    percent.ribo : >= %s",    F$ribo_lo))
message(sprintf("    nFeature_RNA : (%s, %s)", F$nGene_lo, F$nGene_hi))
message(sprintf("    nCount_RNA   : > %s",     F$nUMI_lo))
message(sprintf("  PCA dimensions : %d", NUM_DIM))
message(sprintf("  Chosen CDS     : %s", CHOSEN_CDS))
message(sep)

# ===========================================================================
# 1. Load annotated Seurat object
# ===========================================================================
rds_path <- file.path(p$scratch_dir, "seurat_annotated.rds")
if (!file.exists(rds_path))
  stop("seurat_annotated.rds not found: ", rds_path,
       "\nRun ciri_step02_annotate.R first.", call. = FALSE)

log_info("Loading seurat_annotated.rds ...")
data_seurat <- readRDS(rds_path)
log_info("Loaded: ", ncol(data_seurat), " cells")

pre_QCcounts <- data_seurat[[]] %>%
  group_by(sample_id) %>% summarise(pre_QC_count = n(), .groups = "drop")

# ===========================================================================
# 2. Filter
# ===========================================================================
log_info("Applying filters ...")
filtered_seurat <- subset(
  data_seurat,
  subset = nFeature_RNA  >  F$nGene_lo  &
           nFeature_RNA  <  F$nGene_hi  &
           nCount_RNA    >  F$nUMI_lo   &
           percent.mt    >= F$mito_lo   &
           percent.mt    <  F$mito_hi   &
           percent.ribo  >= F$ribo_lo
)
log_info("Cells after filtering: ", ncol(filtered_seurat))

metadata <- filtered_seurat[[]]
counts   <- SeuratObject::GetAssayData(filtered_seurat, assay = "RNA", layer = "counts")

# Summary table
post_QCcounts <- metadata %>%
  group_by(sample_id) %>% summarise(count = n(), .groups = "drop") %>%
  mutate(
    threshold_nGene = paste0(F$nGene_lo, " – ", F$nGene_hi),
    threshold_nUMI  = paste0("> ",  F$nUMI_lo),
    threshold_mito  = paste0("[",   F$mito_lo, ", ", F$mito_hi, ")%"),
    threshold_ribo  = paste0(">= ", F$ribo_lo, "%")
  ) %>%
  left_join(pre_QCcounts, by = "sample_id") %>%
  mutate(cells_lost = pre_QC_count - count)

write.csv(post_QCcounts,
          paste0(pfx_qc_csv,"_post_QCcounts.csv"), row.names = TRUE)
write.csv(metadata,
          paste0(pfx_qc_csv,"_cell_metadata_filtered_",p$sample,".csv"), row.names = TRUE)
message(""); print(post_QCcounts); message("")

# Post-filter scatter panel
fs1 <- FeatureScatter(filtered_seurat,"nCount_RNA","percent.mt",  group.by="sample_id") +
  geom_vline(xintercept=F$nUMI_lo, linetype="dotted") +
  geom_hline(yintercept=c(F$mito_lo, F$mito_hi), linetype="dotted") +
  scale_colour_manual(values=pal)
fs2 <- FeatureScatter(filtered_seurat,"nCount_RNA","nFeature_RNA", group.by="sample_id") +
  geom_vline(xintercept=F$nUMI_lo, linetype="dotted") +
  geom_hline(yintercept=c(F$nGene_lo, F$nGene_hi), linetype="dotted") +
  scale_colour_manual(values=pal)
fs3 <- FeatureScatter(filtered_seurat,"percent.ribo","percent.mt", group.by="sample_id") +
  geom_vline(xintercept=4, linetype="dotted") +
  geom_hline(yintercept=c(F$mito_lo, F$mito_hi), linetype="dotted") +
  scale_colour_manual(values=pal)
ggsave(paste0(pfx_filter,"_all_QC_after_FILTER.pdf"),
       fs1 + fs2 + fs3, width=20, height=8)

save(metadata, counts, filtered_seurat,
     file = paste0(pfx_qc_R,"_seurat_data_filtered.RData"))

# ===========================================================================
# 3. Monocle3 CDS
# ===========================================================================
log_info("Building Monocle3 CDS ...")
gene_annotation <- data.frame(gene_short_name = rownames(counts),
                               row.names       = rownames(counts))
cds <- new_cell_data_set(counts,
                         cell_metadata = metadata,
                         gene_metadata = gene_annotation)

log_info(sprintf("preprocess_cds(num_dim = %d) ...", NUM_DIM))
cds <- preprocess_cds(cds, num_dim = NUM_DIM)

# ★ Knee plot — inspect to confirm num_dim is appropriate
knee_file <- paste0(pfx_umap,"_variance_knee_plot_dim", NUM_DIM, ".png")
png(knee_file, width=1200, height=800, res=100)
plot_pc_variance_explained(cds)
dev.off()
message("")
message("  ★ INSPECT THE KNEE PLOT: ", knee_file)
message("    If the curve hasn't flattened by dim ", NUM_DIM, ", re-run with --num_dim <X>")
message("")

set.seed(678686)
cds <- reduce_dimension(cds, reduction_method = "tSNE")
cds <- reduce_dimension(cds, umap.min_dist = 0.2, reduction_method = "UMAP")

# ===========================================================================
# 4. Cluster at 5 resolutions
# ===========================================================================
log_info("Clustering at 5 resolutions ...")
res_cfg <- list(
  cds_1 = list(res = 1e-2,  k = NULL, label = "resolution_1e_2"),
  cds_2 = list(res = 1e-4,  k = NULL, label = "resolution_1e_4"),
  cds_3 = list(res = 2e-4,  k = 15,   label = "resolution_2e_4_K15"),
  cds_4 = list(res = 2e-4,  k = NULL, label = "resolution_2e_4"),
  cds_5 = list(res = 5e-4,  k = NULL, label = "resolution_5e_4"),
  cds_6 = list(res = 1e-5,  k = NULL, label = "resolution_1e_5"),
  cds_7 = list(res = 2.5e-5,k = NULL, label = "resolution_25e_6")
)

cds_list <- list()
for (nm in names(res_cfg)) {
  cfg  <- res_cfg[[nm]]
  args <- list(cds, resolution = cfg$res)
  if (!is.null(cfg$k)) args$k <- cfg$k
  clust <- do.call(cluster_cells, args)
  colData(clust)$monocle_clusters <- as.character(monocle3::clusters(clust))
  cds_list[[nm]] <- clust

  p_umap <- plot_cells(clust, color_cells_by="cluster",
                       group_cells_by="cluster", group_label_size=7, cell_size=0.8)
  ggsave(paste0(pfx_umap,"_",nm,"_monocle_clusters_",cfg$label,".pdf"),
         p_umap + facet_grid(~sample_id), width=20, height=4)
  ggsave(paste0(pfx_umap,"_",nm,"_monocle_clusters_",cfg$label,"_single.pdf"),
         p_umap, width=4, height=4)
  log_info("  ", nm, " (res=", cfg$res,
           if(!is.null(cfg$k)) paste0(", k=",cfg$k) else "",
           ") — ", length(unique(monocle3::clusters(clust))), " clusters")
}

cds_1 <- cds_list$cds_1
cds_2 <- cds_list$cds_2
cds_3 <- cds_list$cds_3
cds_4 <- cds_list$cds_4
cds_5 <- cds_list$cds_5
cds_6 <- cds_list$cds_6
cds_7 <- cds_list$cds_7

save(metadata, counts, cds, cds_1, cds_2, cds_3, cds_4, cds_5, cds_6, cds_7,
     file = paste0(pfx_mon_R,"_monocle_cds_clustered.RData"))

# ===========================================================================
# 5. Finalise chosen CDS
# ===========================================================================
log_info("Finalising chosen CDS: ", CHOSEN_CDS)
cds_final <- cds_list[[CHOSEN_CDS]]
cds_final <- learn_graph(cds_final)
cds       <- cds_final

# Add UMAP coords to colData for ggplot outside monocle
pbuild <- plot_cells(cds, color_cells_by="monocle_clusters",
                     label_cell_groups=TRUE, show_trajectory_graph=FALSE,
                     label_leaves=FALSE, label_branch_points=FALSE,
                     graph_label_size=3, group_label_size=4, cell_size=1)
colData(cds)$UMAP_1_monocle <- pbuild$data$data_dim_1
colData(cds)$UMAP_2_monocle <- pbuild$data$data_dim_2

# Final UMAP
p_final <- plot_cells(cds, color_cells_by="cluster", group_cells_by="cluster",
                      group_label_size=7, cell_size=0.8, show_trajectory_graph=FALSE)
ggsave(paste0(pfx_umap,"_UMAP.pdf"),
       p_final + facet_grid(~sample_id), width=10, height=5)
ggsave(paste0(pfx_umap,"_UMAP_single.pdf"), p_final, width=5, height=5)

cell_metadata <- as.data.frame(cds@colData)
save(cell_metadata, counts, cds,
     file = paste0(pfx_mon_R,"_monocle_cds_clustered_FINAL.RData"))
saveRDS(cds, file.path(dirs$mon_scr, "cds_final.rds"))

log_info("Step 02b complete.")
log_info("QC output      : ", qc_root)
log_info("Monocle output : ", mon_root)
message("")
message("═══════════════════════════════════════════════════════════════════")
message("  WHAT TO DO NEXT")
message("═══════════════════════════════════════════════════════════════════")
message("")
message("  1. Open the knee plot:")
message("     ", knee_file)
message("     Curve should flatten before dim ", NUM_DIM, ".")
message("     If it keeps dropping, re-run with a higher --num_dim, e.g.:")
message("")
message("     Rscript ciri_step02_filter.R \\")
message("       --sample    ", p$sample, " \\")
message("       --num_dim   50 \\")
message("       --mito_hi   ", F$mito_hi, "  --mito_lo ", F$mito_lo,
        "  --ribo_lo ", F$ribo_lo, " \\")
message("       --nGene_lo  ", F$nGene_lo, "  --nGene_hi ", F$nGene_hi,
        "  --nUMI_lo ", F$nUMI_lo)
message("")
message("  2. Compare the 7 UMAP resolutions in:")
message("     ", dirs$mon_umap)
message("     _cds_6_ res=1e-5    = fewest clusters  (very coarse)")
message("     _cds_7_ res=2.5e-5  │")
message("     _cds_2_ res=1e-4    │")
message("     _cds_4_ res=2e-4    │")
message("     _cds_3_ res=2e-4 k15│")
message("     _cds_5_ res=5e-4    │")
message("     _cds_1_ res=1e-2    = most clusters   (very fine)")
message("")
message("  3. Once you have chosen the best resolution, run again with")
message("     --chosen_cds, e.g.:")
message("")
message("     Rscript ciri_step02_filter.R \\")
message("       --sample      ", p$sample, " \\")
message("       --chosen_cds  cds_3 \\")
message("       --num_dim     ", NUM_DIM, " \\")
message("       --mito_hi     ", F$mito_hi, "  --mito_lo ", F$mito_lo,
        "  --ribo_lo ", F$ribo_lo, " \\")
message("       --nGene_lo    ", F$nGene_lo, "  --nGene_hi ", F$nGene_hi,
        "  --nUMI_lo ", F$nUMI_lo)
message("")
message("     → saves cds_final.rds to: ", file.path(dirs$mon_scr, "cds_final.rds"))
message("═══════════════════════════════════════════════════════════════════")
