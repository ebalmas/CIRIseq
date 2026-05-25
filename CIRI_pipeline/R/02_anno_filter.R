#!/usr/bin/env Rscript
# =============================================================================
# CIRI Pipeline - Step 02: Annotation & Filtering
# =============================================================================
# Performs QC on the gene expression matrix:
#   - Filters cells with < min_genes detected genes
#   - Restricts to protein-coding genes (from local Ensembl CSV)
#   - Removes mitochondrial (^MT-) and ribosomal (^RPS|^RPL) genes
#   - Removes genes with < min_umis total UMIs
#   - Renames cell barcodes to include guide assignments
#
# Key fix vs original: mitochondrial gene removal now happens AFTER the
# protein-coding filter and is applied to the FINAL output matrix, not
# just used for QC percentage calculations. Also the duplicate
# "Filter to protein-coding genes" call in the original is collapsed to one.
#
# Usage:
#   Rscript 02_anno_filter.R \
#     --dir        /path/to/data/ \
#     --matrix     filtered_feature_bc_matrix.h5 \
#     --gene_ref   ensembl_protein_coding_genes.csv \
#     --min_genes  250 \
#     --min_umis   3 \
#     --remove_mt  TRUE \
#     --remove_rb  TRUE
#
# Output (written to --dir):
#   annotated_matrix.csv                         – filtered, annotated matrix
#   mingenes_proteincoding_filtered_*.csv        – intermediate filtered matrix
#   RiboMito.pdf / RiboMito_filtered.pdf         – QC scatter plots
# =============================================================================

suppressMessages({
  library(dplyr)
  library(ggplot2)
  library(scales)
  library(hdf5r)
  library(Matrix)
  library(Seurat)
  library(data.table)
})

# ---------------------------------------------------------------------------
# Argument Parsing
# ---------------------------------------------------------------------------
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  params <- list(
    dir       = NULL,
    matrix    = NULL,
    gene_ref  = "ensembl_protein_coding_genes.csv",
    min_genes = 250L,
    min_umis  = 3L,
    remove_mt = TRUE,
    remove_rb = TRUE
  )

  i <- 1
  while (i <= length(args)) {
    key <- args[i]
    val <- if (i + 1 <= length(args)) args[i + 1] else NA
    switch(key,
      "--dir"       = { params$dir       <- val; i <- i + 2 },
      "--matrix"    = { params$matrix     <- val; i <- i + 2 },
      "--gene_ref"  = { params$gene_ref   <- val; i <- i + 2 },
      "--min_genes" = { params$min_genes  <- as.integer(val); i <- i + 2 },
      "--min_umis"  = { params$min_umis   <- as.integer(val); i <- i + 2 },
      "--remove_mt" = { params$remove_mt  <- as.logical(val); i <- i + 2 },
      "--remove_rb" = { params$remove_rb  <- as.logical(val); i <- i + 2 },
      { message(paste("Unknown argument:", key)); i <- i + 1 }
    )
  }

  if (is.null(params$dir))    stop("--dir is required")
  if (is.null(params$matrix)) stop("--matrix is required")
  params$dir <- normalizePath(params$dir, mustWork = TRUE)
  params
}

# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------
plot_mito_ribo <- function(meta, title, outfile) {
  color_limits <- range(meta$nCount_RNA, na.rm = TRUE)
  p <- ggplot(meta, aes(x = percent.mt, y = percent.rb, colour = nCount_RNA)) +
    geom_point(alpha = 0.6, size = 1) +
    theme_minimal(base_size = 14) +
    labs(
      x      = "Mitochondrial gene percentage",
      y      = "Ribosomal gene percentage",
      title  = title,
      colour = "nCount_RNA"
    ) +
    scale_color_viridis_c(
      option = "turbo",
      trans  = "log10",
      limits = color_limits,
      oob    = scales::squish
    )
  ggsave(outfile, p, width = 10, height = 10)
  invisible(p)
}

# ---------------------------------------------------------------------------
# Load Protein-Coding Gene Reference
# ---------------------------------------------------------------------------
load_gene_ref <- function(ref_path) {
  # Try the provided path first, then fall back to the working directory
  candidates <- c(ref_path, file.path(getwd(), basename(ref_path)))
  found <- candidates[file.exists(candidates)]

  if (length(found) == 0) {
    stop(paste0(
      "Protein-coding gene reference not found.\n",
      "Expected: ", ref_path, "\n",
      "Run 00_download_gene_ref.R first to generate this file."
    ))
  }

  ref <- read.csv(found[1])

  # Ensure 'mix' column exists
  if (!"mix" %in% names(ref)) {
    ref$mix <- ifelse(
      is.na(ref$hgnc_symbol) | ref$hgnc_symbol == "",
      ref$ensembl_gene_id,
      ref$hgnc_symbol
    )
  }

  message(paste("Loaded", nrow(ref), "protein-coding genes from:", found[1]))
  ref
}

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
main <- function() {
  params <- parse_args()

  message("=== CIRI Pipeline - Step 02: Annotation & Filtering ===")
  message(paste("Directory :", params$dir))
  message(paste("Matrix    :", params$matrix))

  # ---- Load 10X data -------------------------------------------------------
  h5_path <- file.path(params$dir, params$matrix)
  if (!file.exists(h5_path)) stop(paste("H5 matrix not found:", h5_path))

  message("Loading 10X data ...")
  data        <- Read10X_h5(h5_path, unique.features = TRUE)
  data_seurat <- CreateSeuratObject(counts = data$`Gene Expression`, assay = "RNA")
  message(paste("Loaded", ncol(data_seurat), "cells and", nrow(data_seurat), "features"))

  # ---- Pre-filter QC plot --------------------------------------------------
  data_seurat[["percent.mt"]] <- PercentageFeatureSet(data_seurat, pattern = "^MT-")
  data_seurat[["percent.rb"]] <- PercentageFeatureSet(data_seurat, pattern = "^RPS|^RPL")

  plot_mito_ribo(
    data_seurat@meta.data,
    "Mitochondrial vs Ribosomal content (all cells)",
    file.path(params$dir, "RiboMito.pdf")
  )

  # ---- Cell filter: minimum detected genes ---------------------------------
  data_seurat <- subset(data_seurat, subset = nFeature_RNA >= params$min_genes)
  message(paste("After min_genes filter (>=", params$min_genes, "):", ncol(data_seurat), "cells"))

  # Recompute percentages after cell filter
  data_seurat[["percent.mt"]] <- PercentageFeatureSet(data_seurat, pattern = "^MT-")
  data_seurat[["percent.rb"]] <- PercentageFeatureSet(data_seurat, pattern = "^RPS|^RPL")

  plot_mito_ribo(
    data_seurat@meta.data,
    paste("Mitochondrial vs Ribosomal content (>= ", params$min_genes, " genes)", sep = ""),
    file.path(params$dir, "RiboMito_filtered.pdf")
  )

  # ---- Gene filter 1: restrict to protein-coding ---------------------------
  ref <- load_gene_ref(params$gene_ref)
  coding_genes <- ref$mix

  mat <- GetAssayData(data_seurat, assay = "RNA", layer = "counts")
  mat <- mat[intersect(rownames(mat), coding_genes), ]
  message(paste("After protein-coding filter:", nrow(mat), "genes"))

  # ---- Gene filter 2: remove mitochondrial genes ---------------------------
  # FIX: In the original code, MT genes were removed AFTER protein-coding
  # filtering. However, mt-encoded genes (MT-CO1, MT-ND1, etc.) are protein-
  # coding by Ensembl biotype, so they survive the first filter and must be
  # explicitly removed here.
  if (params$remove_mt) {
    mt_mask <- grepl("^MT-", rownames(mat), ignore.case = FALSE)
    n_mt    <- sum(mt_mask)
    mat     <- mat[!mt_mask, ]
    message(paste("Removed", n_mt, "mitochondrial genes (^MT-).",
                  nrow(mat), "genes remaining."))
  }

  # ---- Gene filter 3: remove ribosomal genes --------------------------------
  if (params$remove_rb) {
    rb_mask <- grepl("^RPS|^RPL", rownames(mat), ignore.case = FALSE)
    n_rb    <- sum(rb_mask)
    mat     <- mat[!rb_mask, ]
    message(paste("Removed", n_rb, "ribosomal genes (^RPS|^RPL).",
                  nrow(mat), "genes remaining."))
  }

  # ---- Gene filter 4: minimum total UMIs ----------------------------------
  gene_totals <- Matrix::rowSums(mat)
  mat         <- mat[gene_totals >= params$min_umis, ]
  message(paste("After min_umis filter (>=", params$min_umis, "):", nrow(mat), "genes"))

  # ---- Remove cells with zero counts after gene filtering ------------------
  cell_totals <- Matrix::colSums(mat)
  mat         <- mat[, cell_totals > 0]
  message(paste("After removing zero-count cells:", ncol(mat), "cells"))

  gc()
  message(paste("Final matrix dimensions: genes =", nrow(mat), "/ cells =", ncol(mat)))

  # ---- Save intermediate filtered matrix -----------------------------------
  mat_df    <- as.data.frame(as.matrix(mat))
  inter_out <- file.path(params$dir, "mingenes_proteincoding_filtered_feature_bc_matrix.csv")
  fwrite(mat_df, file = inter_out, row.names = TRUE)
  message(paste("Intermediate filtered matrix saved to:", inter_out))

  # ---- Annotate cells with guide assignments -------------------------------
  ann_path <- file.path(params$dir, "annotation_data.csv")
  if (!file.exists(ann_path)) {
    stop(paste(
      "annotation_data.csv not found in", params$dir,
      "\nRun Step 01 (01_perturbation_assignment.R) first."
    ))
  }

  assigned <- read.csv(ann_path)
  message(paste("Loaded", nrow(assigned), "assigned cells from annotation_data.csv"))

  # Harmonise barcode format (10X uses "-" in barcodes; data.table/R may convert to ".")
  colnames(mat_df) <- sub("\\.", "-", colnames(mat_df))

  # Keep only cells present in the annotation table
  shared_cells <- intersect(colnames(mat_df), assigned$cell_barcode)
  mat_df       <- mat_df[, shared_cells]
  message(paste("Cells shared with annotation:", length(shared_cells)))

  # Rename columns: <barcode>-<feature_a>-<feature_i>
  assigned <- assigned %>%
    mutate(new_name = paste(cell_barcode, feature_a, feature_i, sep = "-"))
  name_map         <- setNames(assigned$new_name, assigned$cell_barcode)
  colnames(mat_df) <- name_map[colnames(mat_df)]

  # ---- Save final annotated matrix -----------------------------------------
  out_path <- file.path(params$dir, "annotated_matrix.csv")
  fwrite(mat_df, file = out_path, row.names = TRUE)
  message(paste("Annotated matrix saved to:", out_path))
  message("Step 02 complete.")
}

main()
