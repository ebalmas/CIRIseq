#!/usr/bin/env Rscript
# =============================================================================
# ciri_step02_annotate.R
# =============================================================================
# Builds an annotated Seurat object and produces QC plots.
# The threshold values below are VISUAL GUIDES ONLY — dotted lines on the
# plots. No cells are removed here. Filtering happens in ciri_step02_filter.R
# after you decide the real cuts from these plots.
#
# INPUTS
# ------
#   scratch/seurat_annotated.rds  (does not exist yet — this script creates it)
#   scratch/filtered_feature_bc_matrix.h5
#   scratch/aggregation.csv
#   scratch/annotation_data.csv           from Step 01
#   scratch/scratch_protein_coding_genes.RData   (optional)
#
# OUTPUT STRUCTURE
# ----------------
#   Output/QC/<YYMMDD>/
#     csv/        <date>_pre_QCcounts.csv
#                 <date>_metadata_annotated.csv
#     ribomito/   <date>_feature_QC.pdf
#                 <date>_UMI_QC.pdf
#                 <date>_gene_det_QC.pdf
#                 <date>_gene_QC.pdf
#                 <date>_gene_UMI_QC.pdf
#                 <date>_percent_MT_RIBO_QC.pdf
#                 <date>_log10GenesPerUMI_QC.pdf
#                 <date>_all_QC_final.pdf
#     R_objects/  <date>_seurat_annotated.RData
#     to_scratch/ seurat_annotated.rds    ← used by ciri_step02_filter.R
#
# USAGE (terminal)
# ----------------
#   Rscript ciri_step02_annotate.R --data_dir /path/to/CIRI_analysis --sample AB014_AB016
#
# USAGE (RStudio — edit INTERACTIVE_PARAMS then source())
# -------------------------------------------------------
#   source("ciri_step02_annotate.R")
# =============================================================================

INTERACTIVE_PARAMS <- list(
  data_dir             = "/home/jovyan/work/shared/CIRI_analysis",
  matrix               = "scratch/filtered_feature_bc_matrix.h5",
  aggr_csv             = "scratch/aggregation.csv",
  scratch_dir          = "scratch",
  output_root          = "Output",
  sample               = "AB014_AB016",
  protein_coding_rdata = "scratch/scratch_protein_coding_genes.RData"  # set "" to skip
)

# ---------------------------------------------------------------------------
# SUGGESTED thresholds — shown as dotted lines on QC plots ONLY.
# Nothing is filtered here. Adjust these to move the dotted lines
# until they visually separate signal from noise, then copy the
# values you like into ciri_step02_filter.R.
# ---------------------------------------------------------------------------
SUGGEST <- list(
  mito_lo  = 1,     # lower mito bound (dotted line)
  mito_hi  = 15,    # upper mito bound (dotted line)
  ribo_lo  = 3,     # lower ribo bound (dotted line)
  nGene_lo = 300,   # min genes per cell (dotted line)
  nGene_hi = 7000,  # max genes per cell (dotted line)
  nUMI_lo  = 100    # min UMIs per cell  (dotted line)
)

# Colour palettes
okabe_pal <- c("#E69F00","#56B4E9","#009E73","#f0E442",
               "#0072B2","#D55E00","#CC79A7","#000000")

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
.parse_args <- function(defaults) {
  raw <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% raw || "-h" %in% raw) {
    cat("Rscript ciri_step02_annotate.R\n",
        "  --data_dir <path>             base directory\n",
        "  --matrix   <rel/abs path>     H5 filename\n",
        "  --aggr_csv <rel/abs path>     aggregation.csv\n",
        "  --scratch_dir <path>          scratch folder\n",
        "  --output_root <path>          Output root\n",
        "  --sample   <name>             experiment name\n",
        "  --protein_coding_rdata <path> optional RData with protein_coding_genes\n",
        "  --mito_lo/hi --ribo_lo        suggested threshold lines on plots\n",
        "  --nGene_lo/hi --nUMI_lo\n")
    quit(status = 0)
  }
  p <- defaults; i <- 1
  while (i <= length(raw)) {
    switch(raw[i],
      "--data_dir"             = { p$data_dir             <- raw[i+1]; i <- i+2 },
      "--matrix"               = { p$matrix                <- raw[i+1]; i <- i+2 },
      "--aggr_csv"             = { p$aggr_csv              <- raw[i+1]; i <- i+2 },
      "--scratch_dir"          = { p$scratch_dir           <- raw[i+1]; i <- i+2 },
      "--output_root"          = { p$output_root           <- raw[i+1]; i <- i+2 },
      "--sample"               = { p$sample                <- raw[i+1]; i <- i+2 },
      "--protein_coding_rdata" = { p$protein_coding_rdata  <- raw[i+1]; i <- i+2 },
      "--mito_lo"              = { SUGGEST$mito_lo  <<- as.numeric(raw[i+1]); i <- i+2 },
      "--mito_hi"              = { SUGGEST$mito_hi  <<- as.numeric(raw[i+1]); i <- i+2 },
      "--ribo_lo"              = { SUGGEST$ribo_lo  <<- as.numeric(raw[i+1]); i <- i+2 },
      "--nGene_lo"             = { SUGGEST$nGene_lo <<- as.numeric(raw[i+1]); i <- i+2 },
      "--nGene_hi"             = { SUGGEST$nGene_hi <<- as.numeric(raw[i+1]); i <- i+2 },
      "--nUMI_lo"              = { SUGGEST$nUMI_lo  <<- as.numeric(raw[i+1]); i <- i+2 },
      { i <- i+1 }
    )
  }
  p
}

p <- .parse_args(INTERACTIVE_PARAMS)

# ---------------------------------------------------------------------------
# Libraries
# ---------------------------------------------------------------------------
suppressMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
})

theme_set(theme_bw(12) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 15, face = "bold", margin = margin(10,0,10,0)),
        axis.text.x = element_text(angle = 45, hjust = 1)))

log_info <- function(...) message("[INFO]  ", ...)
log_warn <- function(...) message("[WARN]  ", ...)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
data_dir <- normalizePath(p$data_dir, mustWork = TRUE)
rp <- function(x) if (startsWith(x, "/")) x else file.path(data_dir, x)

h5_path    <- rp(p$matrix)
aggr_path  <- rp(p$aggr_csv)
ann_path   <- file.path(p$scratch_dir, "annotation_data.csv")

for (f in c(h5_path, aggr_path, ann_path)) {
  if (!file.exists(f)) stop("File not found: ", f, call. = FALSE)
}

# Output dirs
d        <- format(Sys.Date(), "%y%m%d")
qc_root  <- file.path(p$output_root, "QC", d)
dirs <- list(
  csv      = file.path(qc_root, "csv"),
  ribomito = file.path(qc_root, "ribomito"),
  R_obj    = file.path(qc_root, "R_objects"),
  scratch  = file.path(qc_root, "to_scratch")
)
for (dr in dirs) dir.create(dr, recursive = TRUE, showWarnings = FALSE)
pfx_csv <- file.path(dirs$csv,      d)
pfx_rm  <- file.path(dirs$ribomito, d)
pfx_R   <- file.path(dirs$R_obj,    d)

# Banner
sep <- strrep("=", 72)
message(sep)
message("  CIRI — Step 02a: Annotate & QC plots  [NO FILTERING]")
message(sprintf("  %s", Sys.time()))
message(sep)
message("  INPUTS:")
message("    <  ", h5_path)
message("    <  ", aggr_path)
message("    <  ", ann_path)
message("  SUGGESTED dotted lines (not filters):")
message(sprintf("    mito: [%s, %s)   ribo >= %s   nGene: [%s, %s)   nUMI >= %s",
                SUGGEST$mito_lo, SUGGEST$mito_hi, SUGGEST$ribo_lo,
                SUGGEST$nGene_lo, SUGGEST$nGene_hi, SUGGEST$nUMI_lo))
message("  OUTPUTS:")
message("    >  ", qc_root)
message("    >  ", file.path(dirs$scratch, "seurat_annotated.rds"))
message(sep)

# ===========================================================================
# 1. aggregation.csv → sample map (row number = barcode suffix)
# ===========================================================================
log_info("Reading aggregation.csv ...")
aggr_info <- read.csv(aggr_path, stringsAsFactors = FALSE)
aggr_info$orig.ident <- as.character(seq_len(nrow(aggr_info)))

sample_col <- intersect(c("sample_id","library_id","sample"), names(aggr_info))[1]
if (is.na(sample_col))
  stop("aggregation.csv needs 'sample_id' or 'library_id' column. Found: ",
       paste(names(aggr_info), collapse=", "), call. = FALSE)

log_info("Samples:")
for (i in seq_len(nrow(aggr_info)))
  log_info("  suffix -", i, " -> ", aggr_info[[sample_col]][i])

# ===========================================================================
# 2. Load H5 matrix
# ===========================================================================
log_info("Loading H5 matrix ...")
data   <- Read10X_h5(h5_path, unique.features = TRUE)
counts <- if (is.list(data)) {
  log_info("Multiple modalities — using 'Gene Expression'")
  data[["Gene Expression"]]
} else data
log_info("Loaded: ", ncol(counts), " cells, ", nrow(counts), " genes")

# Optional protein-coding filter
pc_rdata <- rp(p$protein_coding_rdata)
if (nchar(p$protein_coding_rdata) > 0 && file.exists(pc_rdata)) {
  log_info("Applying protein-coding gene filter ...")
  env <- new.env()
  load(pc_rdata, envir = env)
  pc_genes <- get(ls(env)[1], envir = env)
  counts   <- counts[rownames(counts) %in% pc_genes, ]
  log_info("  -> ", nrow(counts), " protein-coding genes retained")
}

# ===========================================================================
# 3. Seurat object + aggregation metadata
# ===========================================================================
log_info("Creating Seurat object ...")
data_seurat <- CreateSeuratObject(counts = counts, assay = "RNA")

meta_original          <- data_seurat[[]]
meta_original$libID    <- rownames(meta_original)
meta_original <- meta_original %>%
  mutate(
    orig.ident = sapply(strsplit(libID, "-"), `[[`, 2),
    lib_ID     = sapply(strsplit(libID, "-"), `[[`, 1)
  ) %>%
  right_join(aggr_info, by = "orig.ident")
rownames(meta_original) <- meta_original$libID
data_seurat <- AddMetaData(data_seurat, meta_original)

# ===========================================================================
# 4. Join CIRI guide assignments
# ===========================================================================
log_info("Joining guide assignments from Step 01 ...")
ann <- read.csv(ann_path, stringsAsFactors = FALSE, na.strings = "")

is_a <- !is.na(ann$feature_a) & ann$feature_a != "NA"
is_i <- !is.na(ann$feature_i) & ann$feature_i != "NA"
ann$gene_a    <- ifelse(is_a, sapply(strsplit(ann$feature_a,"_"),`[[`,1), "unassigned")
ann$gene_i    <- ifelse(is_i, sapply(strsplit(ann$feature_i,"_"),`[[`,1), "unassigned")
ann$gene_comb <- paste(ann$gene_a, ann$gene_i, sep="-")
ann$cell_class <- dplyr::case_when(
  is_a & is_i  ~ "CIRI",
  is_a & !is_i ~ "CRISPRa_only",
  !is_a & is_i ~ "CRISPRi_only",
  TRUE         ~ "unassigned"
)
ann$feature_a[!is_a] <- "unassigned"
ann$feature_i[!is_i] <- "unassigned"

ann_join <- ann[, c("cell_barcode","feature_a","feature_i",
                    "gene_a","gene_i","gene_comb","cell_class")]
meta2 <- data_seurat[[]]
meta2$cell_barcode <- rownames(meta2)
meta2 <- left_join(meta2, ann_join, by = "cell_barcode")
for (col in c("feature_a","feature_i","gene_a","gene_i","gene_comb","cell_class"))
  meta2[[col]][is.na(meta2[[col]])] <- "unassigned"
rownames(meta2) <- meta2$cell_barcode
data_seurat <- AddMetaData(data_seurat,
  meta2[, c("feature_a","feature_i","gene_a","gene_i","gene_comb","cell_class")])

log_info("Guide assignment summary:")
print(table(data_seurat$cell_class))

# ===========================================================================
# 5. QC metrics
# ===========================================================================
log_info("Computing QC metrics ...")
DefaultAssay(data_seurat) <- "RNA"
data_seurat[["percent.mt"]]   <- PercentageFeatureSet(data_seurat, pattern = "^MT-")
data_seurat[["percent.ribo"]] <- PercentageFeatureSet(data_seurat, pattern = "RPS")

QC <- data_seurat[[c("nFeature_RNA","percent.mt","nCount_RNA",
                     "percent.ribo","sample_id")]] %>%
  dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA)
QC$info <- rownames(QC)
QC$log10GenesPerUMI <- log10(QC$nGene) / log10(QC$nUMI)

pre_QCcounts <- QC %>% group_by(sample_id) %>% summarise(count = n(), .groups="drop")
write.csv(pre_QCcounts, paste0(pfx_csv,"_pre_QCcounts.csv"), row.names = TRUE)
log_info("Cell counts per sample (pre-filter):"); print(pre_QCcounts)

# Short aliases for readability
S <- SUGGEST

# ===========================================================================
# 6. QC plots  — dotted lines = suggested thresholds, NOT filters
# ===========================================================================
log_info("Saving QC plots ...")

# nGene + nUMI violin
p1 <- ggplot(QC, aes(x=sample_id, y=nGene, fill=sample_id)) +
  geom_point(position=position_jitter(width=0.2), size=0.5, alpha=0.5) +
  geom_violin(scale="width", alpha=0.8) +
  scale_fill_manual(values=okabe_pal) +
  geom_hline(yintercept=c(S$nGene_lo, S$nGene_hi), linetype="dotted") +
  labs(title="nGene per cell  (dotted = suggested cuts)")
p2 <- ggplot(QC, aes(x=sample_id, y=nUMI, fill=sample_id)) +
  geom_point(position=position_jitter(width=0.2), size=0.5, alpha=0.5) +
  geom_violin(scale="width", alpha=0.8) +
  scale_fill_manual(values=okabe_pal) +
  geom_hline(yintercept=S$nUMI_lo, linetype="dotted") +
  labs(title="nUMI per cell  (dotted = suggested cuts)")
ggsave(paste0(pfx_rm,"_feature_QC.pdf"), p1+p2, width=10, height=8)

# UMI density
ggplot(QC, aes(color=sample_id, x=nUMI, fill=sample_id)) +
  geom_density(alpha=0.2) + scale_fill_manual(values=okabe_pal) +
  scale_x_log10() + theme_classic() + ylab("Cell density") +
  geom_vline(xintercept=S$nUMI_lo, linetype="dotted") +
  labs(title="UMI density  (dotted = suggested cut)")
ggsave(paste0(pfx_rm,"_UMI_QC.pdf"), width=10, height=8)

# Gene density
ggplot(QC, aes(color=sample_id, x=nGene, fill=sample_id)) +
  geom_density(alpha=0.2) + scale_fill_manual(values=okabe_pal) +
  theme_classic() + scale_x_log10() +
  geom_vline(xintercept=c(S$nGene_lo, S$nGene_hi), linetype="dotted") +
  labs(title="Gene density  (dotted = suggested cuts)")
ggsave(paste0(pfx_rm,"_gene_det_QC.pdf"), width=10, height=8)

# nGene boxplot
ggplot(QC, aes(x=sample_id, y=log10(nGene), fill=sample_id)) +
  geom_boxplot() + scale_fill_manual(values=okabe_pal) + theme_classic() +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1)) +
  ggtitle("NCells vs NGenes  (dotted = suggested cuts)") +
  geom_hline(yintercept=c(log10(S$nGene_lo), log10(S$nGene_hi)), linetype="dotted")
ggsave(paste0(pfx_rm,"_gene_QC.pdf"), width=10, height=8)

# nUMI vs nGene coloured by mito
ggplot(QC, aes(x=nUMI, y=nGene, color=percent.mt)) +
  geom_point(size=1.5, alpha=0.5) +
  scale_colour_gradient(low="gray90", high="black") +
  stat_smooth(method=lm) + scale_x_log10() + scale_y_log10() + theme_classic() +
  geom_vline(xintercept=S$nUMI_lo, linetype="dotted") +
  geom_hline(yintercept=c(S$nGene_lo, S$nGene_hi), linetype="dotted") +
  labs(title="nUMI vs nGene (colour = %mito; dotted = suggested cuts)")
ggsave(paste0(pfx_rm,"_gene_UMI_QC.pdf"), width=10, height=8)

# nUMI vs nGene by sample
ggplot(QC, aes(x=nUMI, y=nGene, color=sample_id)) +
  geom_point(size=1.5, alpha=0.5) + scale_colour_manual(values=okabe_pal) +
  stat_smooth(method=lm) + scale_x_log10() + scale_y_log10() + theme_classic() +
  geom_vline(xintercept=S$nUMI_lo, linetype="dotted") +
  geom_hline(yintercept=c(S$nGene_lo, S$nGene_hi), linetype="dotted") +
  labs(title="nUMI vs nGene by sample  (dotted = suggested cuts)")
ggsave(paste0(pfx_rm,"_gene_UMI_QC_sample_id.pdf"), width=10, height=8)

# ribo vs mito
ggplot(QC, aes(x=percent.ribo, y=percent.mt, color=sample_id)) +
  geom_point(size=1.5, alpha=0.5) + scale_colour_manual(values=okabe_pal) +
  scale_x_log10() + scale_y_log10() + theme_classic() +
  geom_vline(xintercept=S$ribo_lo,  linetype="dotted") +
  geom_hline(yintercept=c(S$mito_lo, S$mito_hi), linetype="dotted") +
  labs(title="% Ribo vs % Mito  (dotted = suggested cuts)")
ggsave(paste0(pfx_rm,"_percent_MT_RIBO_QC.pdf"), width=10, height=8)

# log10GenesPerUMI
ggplot(QC, aes(x=log10GenesPerUMI, color=sample_id, fill=sample_id)) +
  geom_density(alpha=0.2) + theme_classic() +
  geom_vline(xintercept=0.85, linetype="dotted") +
  labs(title="log10(Genes per UMI)  (dotted = 0.85 suggestion)")
ggsave(paste0(pfx_rm,"_log10GenesPerUMI_QC.pdf"), width=10, height=8)

# FeatureScatter panel
fs1 <- FeatureScatter(data_seurat,"nCount_RNA","percent.mt",  group.by="sample_id") +
  geom_hline(yintercept=c(S$mito_hi, S$mito_lo), linetype="dotted") +
  scale_colour_manual(values=okabe_pal)
fs2 <- FeatureScatter(data_seurat,"nCount_RNA","nFeature_RNA", group.by="sample_id") +
  geom_hline(yintercept=c(S$nGene_lo, S$nGene_hi), linetype="dotted") +
  scale_colour_manual(values=okabe_pal)
fs3 <- FeatureScatter(data_seurat,"percent.ribo","percent.mt", group.by="sample_id") +
  geom_vline(xintercept=S$ribo_lo, linetype="dotted") +
  geom_hline(yintercept=c(S$mito_lo, S$mito_hi), linetype="dotted") +
  scale_colour_manual(values=okabe_pal)
ggsave(paste0(pfx_rm,"_all_QC_final.pdf"), fs1+fs2+fs3, width=20, height=8)

log_info("All QC plots saved to: ", dirs$ribomito)

# ===========================================================================
# 7. Save
# ===========================================================================
write.csv(data_seurat[[]], paste0(pfx_csv,"_metadata_annotated.csv"), row.names=TRUE)
save(meta_original, counts, data_seurat, aggr_info,
     file = paste0(pfx_R,"_seurat_annotated.RData"))
saveRDS(data_seurat, file.path(dirs$scratch, "seurat_annotated.rds"))

log_info("Step 02a complete. Output: ", qc_root)
message("")
message("=========================================================")
message("  NEXT: inspect QC plots in ", dirs$ribomito)
message("  Then open ciri_step02_filter.R, set your real thresholds")
message("  under FILTER_PARAMS, and run.")
message("=========================================================")
