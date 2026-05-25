#!/usr/bin/env Rscript
# =============================================================================
# CIRI Pipeline — Step 02: Annotation & Filtering
# =============================================================================
# INPUT  (auto-resolved from Step 01 to_scratch/, or override with --scratch):
#   annotation_data.csv
# INPUT  (provide directly):
#   --data_dir   folder with the H5 matrix
#   --gene_ref   ensembl_protein_coding_genes.csv (from Step 00 to_scratch/)
#
# OUTPUT STRUCTURE:
#   Output/
#     <YYMMDD>_step02_filter_<sample>/
#       plots/      RiboMito.pdf, RiboMito_filtered.pdf
#       stats/      filter_summary.txt
#       R_objects/  (none)
#       to_scratch/ annotated_matrix.csv   ← consumed by Step 03
#
# Usage:
#   Rscript steps/step_02_anno_filter.R \
#     --data_dir    /path/to/data \
#     --matrix      filtered_feature_bc_matrix.h5 \
#     --sample      AB011
#   Rscript steps/step_02_anno_filter.R --help
# =============================================================================

STEPS_DIR <- dirname(sys.frame(1)$ofile)
source(file.path(STEPS_DIR, "utils.R"))

SPEC <- list(
  data_dir    = list(flag="--data_dir",    required=TRUE,  type="character",
                     help="Directory containing the H5 matrix"),
  matrix      = list(flag="--matrix",      required=TRUE,  type="character",
                     help="H5 matrix filename inside --data_dir"),
  output_root = list(flag="--output_root", required=FALSE, type="character",
                     default="Output",    help="Top-level output directory"),
  sample      = list(flag="--sample",      required=FALSE, type="character",
                     default="CIRI",      help="Sample / experiment name"),
  scratch     = list(flag="--scratch",     required=FALSE, type="character",
                     default=NULL,
                     help="Path to Step 01 to_scratch/ folder (auto-resolved if omitted)"),
  gene_ref    = list(flag="--gene_ref",    required=FALSE, type="character",
                     default=NULL,
                     help="Path to ensembl_protein_coding_genes.csv (auto-resolved from Step 00 to_scratch/ if omitted)"),
  min_genes   = list(flag="--min_genes",   required=FALSE, type="integer",
                     default=250L,        help="Min detected genes per cell"),
  min_umis    = list(flag="--min_umis",    required=FALSE, type="integer",
                     default=3L,          help="Min total UMIs per gene"),
  remove_mt   = list(flag="--remove_mt",   required=FALSE, type="logical",
                     default=TRUE,        help="Remove mitochondrial genes (^MT-)"),
  remove_rb   = list(flag="--remove_rb",   required=FALSE, type="logical",
                     default=TRUE,        help="Remove ribosomal genes (^RPS|^RPL)")
)

main <- function() {
  suppressMessages({
    library(ggplot2); library(scales); library(hdf5r)
    library(Matrix); library(Seurat); library(data.table); library(dplyr)
  })

  p   <- parse_args(SPEC)
  out <- make_out_dirs(p$output_root, "step02_filter", p$sample)

  # Resolve annotation from Step 01 to_scratch/
  scratch01 <- if (!is.null(p$scratch)) p$scratch else
    find_scratch(p$output_root, "step01_assignment", p$sample)
  ann_path  <- file.path(scratch01, "annotation_data.csv")

  # Resolve gene reference from Step 00 to_scratch/
  gene_ref_path <- if (!is.null(p$gene_ref)) {
    p$gene_ref
  } else {
    scratch00 <- tryCatch(find_scratch(p$output_root, "step00_geneRef", p$sample),
                          error = function(e) NULL)
    if (!is.null(scratch00)) file.path(scratch00, "ensembl_protein_coding_genes.csv") else NULL
  }

  step_banner("02", "Annotation & Filtering",
    inputs  = c(file.path(p$data_dir, p$matrix), ann_path,
                if (!is.null(gene_ref_path)) gene_ref_path else "ensembl_protein_coding_genes.csv [not found]"),
    outputs = c(out$plots, out$stats,
                file.path(out$to_scratch, "annotated_matrix.csv"))
  )

  p$data_dir <- normalizePath(p$data_dir, mustWork = TRUE)
  assert_file(file.path(p$data_dir, p$matrix))
  assert_file(ann_path, hint = "Run step_01_perturbation_assignment.R first.")
  if (is.null(gene_ref_path) || !file.exists(gene_ref_path))
    stop("Ensembl gene reference not found.\nRun step_00_download_gene_ref.R first, or pass --gene_ref.")

  # ---- Load 10X data -----------------------------------------------------
  log_info("Loading 10X data ...")
  data <- Read10X_h5(file.path(p$data_dir, p$matrix), unique.features = TRUE)
  seu  <- CreateSeuratObject(counts = data$`Gene Expression`, assay = "RNA")
  log_info("Loaded: ", ncol(seu), " cells, ", nrow(seu), " features")

  mito_ribo_plot <- function(seu, title, outfile) {
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
    seu[["percent.rb"]] <- PercentageFeatureSet(seu, pattern = "^RPS|^RPL")
    lims <- range(seu@meta.data$nCount_RNA, na.rm = TRUE)
    p <- ggplot(seu@meta.data, aes(x=percent.mt, y=percent.rb, colour=nCount_RNA)) +
      geom_point(alpha=0.6, size=1) + theme_minimal(base_size=14) +
      labs(x="% Mitochondrial", y="% Ribosomal", title=title, colour="nCount_RNA") +
      scale_color_viridis_c(option="turbo", trans="log10", limits=lims, oob=scales::squish)
    ggsave(outfile, p, width=10, height=10)
  }

  mito_ribo_plot(seu, "Mito vs Ribo — all cells",
                 file.path(out$plots, "RiboMito_pre_filter.pdf"))

  # ---- Cell filter -------------------------------------------------------
  n_before <- ncol(seu)
  seu      <- subset(seu, subset = nFeature_RNA >= p$min_genes)
  log_info("Cell filter (>= ", p$min_genes, " genes): ", n_before, " → ", ncol(seu))
  mito_ribo_plot(seu, paste0("Mito vs Ribo — >= ", p$min_genes, " genes"),
                 file.path(out$plots, "RiboMito_post_cell_filter.pdf"))

  # ---- Gene filter 1: protein-coding -------------------------------------
  ref <- read.csv(gene_ref_path)
  if (!"mix" %in% names(ref))
    ref$mix <- ifelse(is.na(ref$hgnc_symbol) | ref$hgnc_symbol == "",
                      ref$ensembl_gene_id, ref$hgnc_symbol)
  mat    <- GetAssayData(seu, assay="RNA", layer="counts")
  n_all  <- nrow(mat)
  mat    <- mat[intersect(rownames(mat), ref$mix), ]
  log_info("Protein-coding filter: ", n_all, " → ", nrow(mat), " genes")

  # ---- Gene filter 2: remove mitochondrial (FIX) -------------------------
  n_mt <- 0
  if (p$remove_mt) {
    mt  <- grepl("^MT-", rownames(mat))
    n_mt <- sum(mt); mat <- mat[!mt, ]
    log_info("Removed ", n_mt, " mitochondrial genes → ", nrow(mat))
  }

  # ---- Gene filter 3: remove ribosomal -----------------------------------
  n_rb <- 0
  if (p$remove_rb) {
    rb  <- grepl("^RPS|^RPL", rownames(mat))
    n_rb <- sum(rb); mat <- mat[!rb, ]
    log_info("Removed ", n_rb, " ribosomal genes → ", nrow(mat))
  }

  # ---- Gene filter 4: min UMIs -------------------------------------------
  gene_umi <- Matrix::rowSums(mat)
  n_low    <- sum(gene_umi < p$min_umis)
  mat      <- mat[gene_umi >= p$min_umis, ]
  log_info("Low-UMI gene filter (< ", p$min_umis, "): removed ", n_low, " → ", nrow(mat))

  cell_umi <- Matrix::colSums(mat)
  mat      <- mat[, cell_umi > 0]
  log_info("Final matrix: ", nrow(mat), " genes × ", ncol(mat), " cells")

  # ---- Annotate barcodes -------------------------------------------------
  assigned <- read.csv(ann_path)
  mat_df   <- as.data.frame(as.matrix(mat))
  colnames(mat_df) <- sub("\\.", "-", colnames(mat_df))
  shared   <- intersect(colnames(mat_df), assigned$cell_barcode)
  mat_df   <- mat_df[, shared]
  log_info("Cells shared with annotation: ", length(shared))
  assigned <- dplyr::mutate(assigned,
    new_name = paste(cell_barcode, feature_a, feature_i, sep="-"))
  name_map         <- setNames(assigned$new_name, assigned$cell_barcode)
  colnames(mat_df) <- name_map[colnames(mat_df)]

  # ---- stats/ ------------------------------------------------------------
  stats_txt <- c(
    paste("Sample            :", p$sample),
    paste("Input cells       :", n_before),
    paste("After min_genes   :", ncol(seu)),
    paste("MT genes removed  :", n_mt),
    paste("RB genes removed  :", n_rb),
    paste("Low-UMI genes rm  :", n_low),
    paste("Final genes       :", nrow(mat_df)),
    paste("Final cells       :", ncol(mat_df))
  )
  writeLines(stats_txt, file.path(out$stats, "filter_summary.txt"))
  message(paste(stats_txt, collapse="\n"))

  # ---- to_scratch/ -------------------------------------------------------
  data.table::fwrite(mat_df,
    file = file.path(out$to_scratch, "annotated_matrix.csv"), row.names = TRUE)

  log_info("Step 02 complete.")
  log_info("Output folder: ", out$out_dir)
}

main()
