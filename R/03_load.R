#!/usr/bin/env Rscript
# =============================================================================
# CIRI Pipeline - Step 03: Data Loading & Preprocessing
# =============================================================================
# Initializes a Monocle3 cell_data_set object, performs normalization,
# PCA, and generates the initial UMAP.
#
# Usage:
#   Rscript 03_load.R \
#     --dir         /path/to/data/ \
#     --matrix      annotated_matrix.csv \
#     --resolution  0.5e-4 \
#     --seed        1234597698
#
# Output (written to --dir):
#   processed_cds.RData    – Monocle3 cds object
#   *.pdf                  – UMAP plots coloured by cluster, guide, sample
# =============================================================================

suppressMessages({
  library(monocle3)
  library(tidyverse)
  library(ggplot2)
  library(gtools)
  library(dplyr)
  library(Seurat)
  library(viridis)
  library(ggrepel)
})

# ---------------------------------------------------------------------------
# Argument Parsing
# ---------------------------------------------------------------------------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  params <- list(
    dir        = NULL,
    matrix     = "annotated_matrix.csv",
    resolution = 0.5e-4,
    seed       = 1234597698L
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    val <- if (i + 1 <= length(args)) args[i + 1] else NA
    switch(key,
      "--dir"        = { params$dir        <- val; i <- i + 2 },
      "--matrix"     = { params$matrix      <- val; i <- i + 2 },
      "--resolution" = { params$resolution  <- as.numeric(val); i <- i + 2 },
      "--seed"       = { params$seed        <- as.integer(val); i <- i + 2 },
      { message(paste("Unknown argument:", key)); i <- i + 1 }
    )
  }

  if (is.null(params$dir)) stop("--dir is required")
  params$dir <- normalizePath(params$dir, mustWork = TRUE)
  params
}

# ---------------------------------------------------------------------------
# Build cell metadata from column names
# Format: <barcode>-<sample>-<feature_a>-<feature_i>
# ---------------------------------------------------------------------------
build_cell_meta <- function(col_names) {
  data <- data.frame(nomi = col_names, row.names = col_names, stringsAsFactors = FALSE)

  # Detect if cell names include sample index (4+ dash-separated fields)
  data <- tidyr::separate(
    data, nomi,
    into   = c("cellID", "sample", "guide_a", "guide_i"),
    sep    = "-",
    remove = FALSE,
    fill   = "right"
  )

  data <- data %>%
    mutate(
      guide_a = ifelse(guide_a %in% c("NA", "NA;NA", NA), NA_character_, guide_a),
      guide_i = ifelse(guide_i %in% c("NA", "NA;NA", NA), NA_character_, guide_i),
      comb    = paste(guide_a, guide_i, sep = "-"),
      gene_a  = sapply(strsplit(guide_a, "_"), `[`, 1),
      gene_i  = sapply(strsplit(guide_i, "_"), `[`, 1),
      gene_comb = paste(gene_a, gene_i, sep = "-")
    )

  data
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
main <- function() {
  params <- parse_args()
  set.seed(params$seed)

  message("=== CIRI Pipeline - Step 03: Data Loading & Preprocessing ===")
  message(paste("Directory  :", params$dir))
  message(paste("Matrix     :", params$matrix))
  message(paste("Resolution :", params$resolution))

  # ---- Load expression matrix ----------------------------------------------
  mat_path <- file.path(params$dir, params$matrix)
  if (!file.exists(mat_path)) stop(paste("Matrix file not found:", mat_path))

  message("Reading expression matrix ...")
  exp         <- read.csv(mat_path, header = TRUE, check.names = FALSE)
  rownames(exp) <- exp[, 1]
  exp[, 1]    <- NULL
  message(paste("Matrix dimensions: genes =", nrow(exp), "/ cells =", ncol(exp)))

  # ---- Build cell metadata ------------------------------------------------
  cell_meta <- build_cell_meta(colnames(exp))

  # ---- Build Monocle3 CDS -------------------------------------------------
  message("Creating Monocle3 cell_data_set ...")
  gene_meta <- data.frame(
    gene_short_name = rownames(exp),
    row.names       = rownames(exp)
  )

  cds <- new_cell_data_set(
    expression_data = as.matrix(exp),
    cell_metadata   = cell_meta,
    gene_metadata   = gene_meta
  )

  # ---- Preprocess (normalization + PCA) -----------------------------------
  message("Preprocessing (normalization, PCA) ...")
  cds <- preprocess_cds(cds, num_dim = 100)

  # ---- Alignment (optional: removes batch effects between samples) ---------
  # Uncomment and set 'residual_model_formula_str' if you have multiple samples
  # cds <- align_cds(cds, alignment_group = "sample")

  # ---- UMAP ---------------------------------------------------------------
  message("Running UMAP ...")
  cds <- reduce_dimension(cds)

  # ---- Clustering ----------------------------------------------------------
  message(paste("Clustering (resolution:", params$resolution, ") ..."))
  cds <- cluster_cells(cds, resolution = params$resolution)

  # ---- Save plots ----------------------------------------------------------
  message("Saving UMAP plots ...")

  p_cluster <- plot_cells(cds, color_cells_by = "cluster",
                          label_cell_groups = TRUE, label_leaves = FALSE,
                          label_branch_points = FALSE, graph_label_size = 2) +
    ggtitle("UMAP — Cluster")
  ggsave(file.path(params$dir, "umap_clusters.pdf"), p_cluster, width = 10, height = 8)

  p_sample <- plot_cells(cds, color_cells_by = "sample",
                         label_cell_groups = FALSE) +
    ggtitle("UMAP — Sample")
  ggsave(file.path(params$dir, "umap_sample.pdf"), p_sample, width = 10, height = 8)

  if ("gene_a" %in% names(colData(cds))) {
    p_guide_a <- plot_cells(cds, color_cells_by = "gene_a",
                            label_cell_groups = FALSE) +
      ggtitle("UMAP — CRISPRa gene")
    ggsave(file.path(params$dir, "umap_guide_a.pdf"), p_guide_a, width = 10, height = 8)
  }

  if ("gene_i" %in% names(colData(cds))) {
    p_guide_i <- plot_cells(cds, color_cells_by = "gene_i",
                            label_cell_groups = FALSE) +
      ggtitle("UMAP — CRISPRi gene")
    ggsave(file.path(params$dir, "umap_guide_i.pdf"), p_guide_i, width = 10, height = 8)
  }

  # ---- Save CDS object ----------------------------------------------------
  out_rdata <- file.path(params$dir, "processed_cds.RData")
  save(cds, file = out_rdata)
  message(paste("CDS saved to:", out_rdata))
  message("Step 03 complete.")
}

main()
