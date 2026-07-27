# =============================================================================
# CIRI — Step 02a: Annotate & QC   |   Step 02b: Filter & Monocle3
# =============================================================================
# These are the package-function equivalents of the standalone scripts
# ciri_step02_annotate.R and ciri_step02_filter.R.
# The standalone scripts in inst/scripts/ remain available for direct use.
# =============================================================================

# ---- Shared colour palettes (used by both functions) -----------------------
.okabe_pal <- c("#E69F00","#56B4E9","#009E73","#f0E442",
                "#0072B2","#D55E00","#CC79A7","#000000")

.pal15 <- c("#000000","#004949","#009292","#ff6db6","#ffb6db",
            "#490092","#006ddb","#b66dff","#6db6ff","#b6dbff",
            "#920000","#924900","#db6d00","#24ff24","#ffff6d")

# ---- Shared ggplot theme ---------------------------------------------------
.ciri_theme <- function() {
  ggplot2::theme_bw(12) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title       = ggplot2::element_text(size = 15, face = "bold",
                                               margin = ggplot2::margin(10,0,10,0)),
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1)
    )
}

# =============================================================================
# ciri_step02_annotate
# =============================================================================

#' Build an annotated Seurat object and produce QC plots
#'
#' Reads the CellRanger-aggr H5 matrix, joins guide assignments from Step 01
#' and sample information from \code{aggregation.csv}, computes mito/ribo
#' percentages, and saves QC plots with suggested threshold dotted lines.
#' \strong{No cells are removed.} Filtering is done separately by
#' \code{\link{ciri_step02_filter}}.
#'
#' @param data_dir             Base directory containing the H5 and aggregation CSV.
#' @param matrix               H5 filename, relative to \code{data_dir}.
#'                             Default \code{"scratch/filtered_feature_bc_matrix.h5"}.
#' @param aggr_csv             aggregation.csv filename, relative to \code{data_dir}.
#'                             Default \code{"scratch/aggregation.csv"}.
#' @param scratch_dir          Folder containing \code{annotation_data.csv} from
#'                             Step 01. Default \code{"scratch"}.
#' @param output_root          Top-level output directory. Default \code{"Output"}.
#' @param sample               Experiment name used in output folder names.
#' @param protein_coding_rdata Optional path to an \code{.RData} file containing
#'                             a character vector of protein-coding gene symbols.
#'                             If provided, the matrix is subsetted to these genes.
#'                             Set \code{""} to skip.
#' @param suggest_mito_lo      Suggested lower mito threshold shown as a dotted
#'                             line on QC plots. Default \code{1}.
#' @param suggest_mito_hi      Suggested upper mito threshold. Default \code{15}.
#' @param suggest_ribo_lo      Suggested lower ribo threshold. Default \code{3}.
#' @param suggest_nGene_lo     Suggested min genes per cell. Default \code{300}.
#' @param suggest_nGene_hi     Suggested max genes per cell. Default \code{7000}.
#' @param suggest_nUMI_lo      Suggested min UMIs per cell. Default \code{100}.
#'
#' @return Invisibly returns the path to the output folder.
#'
#' @section Output structure:
#' \preformatted{
#' Output/QC/<YYMMDD>/
#'   csv/        <date>_pre_QCcounts.csv
#'               <date>_metadata_annotated.csv
#'   ribomito/   <date>_feature_QC.pdf
#'               <date>_UMI_QC.pdf
#'               <date>_gene_det_QC.pdf
#'               <date>_gene_QC.pdf
#'               <date>_gene_UMI_QC.pdf
#'               <date>_percent_MT_RIBO_QC.pdf
#'               <date>_log10GenesPerUMI_QC.pdf
#'               <date>_all_QC_final.pdf
#'   R_objects/  <date>_seurat_annotated.RData
#'   to_scratch/ seurat_annotated.rds    <- used by ciri_step02_filter()
#' }
#'
#' @examples
#' \dontrun{
#' ciri_step02_annotate(
#'   data_dir  = "/path/to/CIRI_analysis",
#'   sample    = "AB014_AB016",
#'   suggest_mito_hi  = 15,
#'   suggest_nGene_lo = 300
#' )
#' }
#' @export
ciri_step02_annotate <- function(data_dir,
                                  matrix               = "scratch/filtered_feature_bc_matrix.h5",
                                  aggr_csv             = "scratch/aggregation.csv",
                                  scratch_dir          = "scratch",
                                  output_root          = "Output",
                                  sample               = "CIRI",
                                  protein_coding_rdata = "",
                                  suggest_mito_lo      = 1,
                                  suggest_mito_hi      = 15,
                                  suggest_ribo_lo      = 3,
                                  suggest_nGene_lo     = 300,
                                  suggest_nGene_hi     = 7000,
                                  suggest_nUMI_lo      = 100) {

  .check_cran_pkgs("Seurat", "SeuratObject", "patchwork")

  data_dir <- normalizePath(data_dir, mustWork = TRUE)
  rp       <- function(x) if (startsWith(x, "/")) x else file.path(data_dir, x)

  h5_path   <- rp(matrix)
  aggr_path <- rp(aggr_csv)
  ann_path  <- file.path(scratch_dir, "annotation_data.csv")

  for (f in c(h5_path, aggr_path, ann_path))
    assert_file(f)

  # Output dirs
  d       <- format(Sys.Date(), "%y%m%d")
  qc_root <- file.path(output_root, "QC", d)
  dirs    <- list(
    csv      = file.path(qc_root, "csv"),
    ribomito = file.path(qc_root, "ribomito"),
    R_obj    = file.path(qc_root, "R_objects"),
    scratch  = file.path(qc_root, "to_scratch")
  )
  for (dr in dirs) dir.create(dr, recursive = TRUE, showWarnings = FALSE)
  pfx_csv <- file.path(dirs$csv,      d)
  pfx_rm  <- file.path(dirs$ribomito, d)
  pfx_R   <- file.path(dirs$R_obj,    d)

  S <- list(mito_lo = suggest_mito_lo, mito_hi = suggest_mito_hi,
            ribo_lo = suggest_ribo_lo,
            nGene_lo = suggest_nGene_lo, nGene_hi = suggest_nGene_hi,
            nUMI_lo = suggest_nUMI_lo)

  step_banner("02a", "Annotate & QC  [NO FILTERING]",
    inputs  = c(h5_path, aggr_path, ann_path),
    outputs = c(qc_root, file.path(dirs$scratch, "seurat_annotated.rds"))
  )
  message(sprintf("  Suggested dotted lines: mito [%s,%s)  ribo>=%s  nGene [%s,%s)  nUMI>=%s",
                  S$mito_lo, S$mito_hi, S$ribo_lo,
                  S$nGene_lo, S$nGene_hi, S$nUMI_lo))

  # ---- 1. aggregation.csv -------------------------------------------------
  log_info("Reading aggregation.csv ...")
  aggr_info           <- utils::read.csv(aggr_path, stringsAsFactors = FALSE)
  aggr_info$orig.ident <- as.character(seq_len(nrow(aggr_info)))
  sample_col          <- intersect(c("sample_id","library_id","sample"),
                                   names(aggr_info))[1]
  if (is.na(sample_col))
    stop("aggregation.csv needs 'sample_id' or 'library_id' column.\nFound: ",
         paste(names(aggr_info), collapse=", "), call. = FALSE)
  for (i in seq_len(nrow(aggr_info)))
    log_info("  suffix -", i, " -> ", aggr_info[[sample_col]][i])

  # ---- 2. H5 matrix -------------------------------------------------------
  log_info("Loading H5 matrix ...")
  data   <- Seurat::Read10X_h5(h5_path, unique.features = TRUE)
  counts <- if (is.list(data)) {
    log_info("Multiple modalities — using 'Gene Expression'")
    data[["Gene Expression"]]
  } else data
  log_info("Loaded: ", ncol(counts), " cells, ", nrow(counts), " genes")

  if (nchar(protein_coding_rdata) > 0 && file.exists(rp(protein_coding_rdata))) {
    log_info("Applying protein-coding gene filter ...")
    env      <- new.env()
    load(rp(protein_coding_rdata), envir = env)
    pc_genes <- get(ls(env)[1], envir = env)
    counts   <- counts[rownames(counts) %in% pc_genes, ]
    log_info("  -> ", nrow(counts), " protein-coding genes retained")
  }

  # ---- 3. Seurat + aggregation metadata -----------------------------------
  log_info("Creating Seurat object ...")
  data_seurat            <- Seurat::CreateSeuratObject(counts = counts, assay = "RNA")
  meta_original          <- data_seurat[[]]
  meta_original$libID    <- rownames(meta_original)
  meta_original <- meta_original |>
    dplyr::mutate(
      orig.ident = sapply(strsplit(libID, "-"), `[[`, 2),
      lib_ID     = sapply(strsplit(libID, "-"), `[[`, 1)
    ) |>
    dplyr::right_join(aggr_info, by = "orig.ident")
  rownames(meta_original) <- meta_original$libID
  data_seurat <- Seurat::AddMetaData(data_seurat, meta_original)

  # ---- 4. Guide assignments -----------------------------------------------
  log_info("Joining guide assignments from Step 01 ...")
  ann   <- utils::read.csv(ann_path, stringsAsFactors = FALSE, na.strings = "")
  is_a  <- !is.na(ann$feature_a) & ann$feature_a != "NA"
  is_i  <- !is.na(ann$feature_i) & ann$feature_i != "NA"
  ann$gene_a    <- ifelse(is_a, sapply(strsplit(ann$feature_a,"_"),`[[`,1), "unassigned")
  ann$gene_i    <- ifelse(is_i, sapply(strsplit(ann$feature_i,"_"),`[[`,1), "unassigned")
  ann$gene_comb <- paste(ann$gene_a, ann$gene_i, sep = "-")
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
  meta2              <- data_seurat[[]]
  meta2$cell_barcode <- rownames(meta2)
  meta2 <- dplyr::left_join(meta2, ann_join, by = "cell_barcode")
  for (col in c("feature_a","feature_i","gene_a","gene_i","gene_comb","cell_class"))
    meta2[[col]][is.na(meta2[[col]])] <- "unassigned"
  rownames(meta2) <- meta2$cell_barcode
  data_seurat <- Seurat::AddMetaData(
    data_seurat,
    meta2[, c("feature_a","feature_i","gene_a","gene_i","gene_comb","cell_class")]
  )
  log_info("Guide assignment summary:")
  print(table(data_seurat$cell_class))

  # ---- 5. QC metrics ------------------------------------------------------
  log_info("Computing QC metrics ...")
  Seurat::DefaultAssay(data_seurat) <- "RNA"
  data_seurat[["percent.mt"]]   <- Seurat::PercentageFeatureSet(data_seurat, pattern = "^MT-")
  data_seurat[["percent.ribo"]] <- Seurat::PercentageFeatureSet(data_seurat, pattern = "RPS")

  QC <- data_seurat[[c("nFeature_RNA","percent.mt","nCount_RNA",
                        "percent.ribo","sample_id")]] |>
    dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA)
  QC$info              <- rownames(QC)
  QC$log10GenesPerUMI  <- log10(QC$nGene) / log10(QC$nUMI)

  pre_QCcounts <- QC |> dplyr::group_by(sample_id) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop")
  utils::write.csv(pre_QCcounts, paste0(pfx_csv,"_pre_QCcounts.csv"), row.names = TRUE)
  log_info("Cell counts per sample (pre-filter):"); print(pre_QCcounts)

  # ---- 6. QC plots --------------------------------------------------------
  log_info("Saving QC plots ...")

  p1 <- ggplot2::ggplot(QC, ggplot2::aes(x=sample_id, y=nGene, fill=sample_id)) +
    ggplot2::geom_point(position=ggplot2::position_jitter(width=0.2),size=0.5,alpha=0.5) +
    ggplot2::geom_violin(scale="width", alpha=0.8) +
    ggplot2::scale_fill_manual(values=.okabe_pal) +
    ggplot2::geom_hline(yintercept=c(S$nGene_lo,S$nGene_hi), linetype="dotted") +
    ggplot2::labs(title="nGene per cell  (dotted = suggested cuts)") + .ciri_theme()
  p2 <- ggplot2::ggplot(QC, ggplot2::aes(x=sample_id, y=nUMI, fill=sample_id)) +
    ggplot2::geom_point(position=ggplot2::position_jitter(width=0.2),size=0.5,alpha=0.5) +
    ggplot2::geom_violin(scale="width", alpha=0.8) +
    ggplot2::scale_fill_manual(values=.okabe_pal) +
    ggplot2::geom_hline(yintercept=S$nUMI_lo, linetype="dotted") +
    ggplot2::labs(title="nUMI per cell  (dotted = suggested cuts)") + .ciri_theme()
  ggplot2::ggsave(paste0(pfx_rm,"_feature_QC.pdf"), patchwork::wrap_plots(p1, p2), width=10, height=8)

  ggplot2::ggplot(QC, ggplot2::aes(color=sample_id, x=nUMI, fill=sample_id)) +
    ggplot2::geom_density(alpha=0.2) + ggplot2::scale_fill_manual(values=.okabe_pal) +
    ggplot2::scale_x_log10() + ggplot2::theme_classic() + ggplot2::ylab("Cell density") +
    ggplot2::geom_vline(xintercept=S$nUMI_lo, linetype="dotted") +
    ggplot2::labs(title="UMI density  (dotted = suggested cut)")
  ggplot2::ggsave(paste0(pfx_rm,"_UMI_QC.pdf"), width=10, height=8)

  ggplot2::ggplot(QC, ggplot2::aes(color=sample_id, x=nGene, fill=sample_id)) +
    ggplot2::geom_density(alpha=0.2) + ggplot2::scale_fill_manual(values=.okabe_pal) +
    ggplot2::theme_classic() + ggplot2::scale_x_log10() +
    ggplot2::geom_vline(xintercept=c(S$nGene_lo,S$nGene_hi), linetype="dotted") +
    ggplot2::labs(title="Gene density  (dotted = suggested cuts)")
  ggplot2::ggsave(paste0(pfx_rm,"_gene_det_QC.pdf"), width=10, height=8)

  ggplot2::ggplot(QC, ggplot2::aes(x=sample_id, y=log10(nGene), fill=sample_id)) +
    ggplot2::geom_boxplot() + ggplot2::scale_fill_manual(values=.okabe_pal) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=45, vjust=1, hjust=1)) +
    ggplot2::ggtitle("NCells vs NGenes  (dotted = suggested cuts)") +
    ggplot2::geom_hline(yintercept=c(log10(S$nGene_lo),log10(S$nGene_hi)), linetype="dotted")
  ggplot2::ggsave(paste0(pfx_rm,"_gene_QC.pdf"), width=10, height=8)

  ggplot2::ggplot(QC, ggplot2::aes(x=nUMI, y=nGene, color=percent.mt)) +
    ggplot2::geom_point(size=1.5, alpha=0.5) +
    ggplot2::scale_colour_gradient(low="gray90", high="black") +
    ggplot2::stat_smooth(method=stats::lm) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() + ggplot2::theme_classic() +
    ggplot2::geom_vline(xintercept=S$nUMI_lo, linetype="dotted") +
    ggplot2::geom_hline(yintercept=c(S$nGene_lo,S$nGene_hi), linetype="dotted") +
    ggplot2::labs(title="nUMI vs nGene (colour = %mito; dotted = suggested cuts)")
  ggplot2::ggsave(paste0(pfx_rm,"_gene_UMI_QC.pdf"), width=10, height=8)

  ggplot2::ggplot(QC, ggplot2::aes(x=nUMI, y=nGene, color=sample_id)) +
    ggplot2::geom_point(size=1.5, alpha=0.5) +
    ggplot2::scale_colour_manual(values=.okabe_pal) +
    ggplot2::stat_smooth(method=stats::lm) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() + ggplot2::theme_classic() +
    ggplot2::geom_vline(xintercept=S$nUMI_lo, linetype="dotted") +
    ggplot2::geom_hline(yintercept=c(S$nGene_lo,S$nGene_hi), linetype="dotted") +
    ggplot2::labs(title="nUMI vs nGene by sample  (dotted = suggested cuts)")
  ggplot2::ggsave(paste0(pfx_rm,"_gene_UMI_QC_sample_id.pdf"), width=10, height=8)

  ggplot2::ggplot(QC, ggplot2::aes(x=percent.ribo, y=percent.mt, color=sample_id)) +
    ggplot2::geom_point(size=1.5, alpha=0.5) +
    ggplot2::scale_colour_manual(values=.okabe_pal) +
    ggplot2::scale_x_log10() + ggplot2::scale_y_log10() + ggplot2::theme_classic() +
    ggplot2::geom_vline(xintercept=S$ribo_lo, linetype="dotted") +
    ggplot2::geom_hline(yintercept=c(S$mito_lo,S$mito_hi), linetype="dotted") +
    ggplot2::labs(title="% Ribo vs % Mito  (dotted = suggested cuts)")
  ggplot2::ggsave(paste0(pfx_rm,"_percent_MT_RIBO_QC.pdf"), width=10, height=8)

  ggplot2::ggplot(QC, ggplot2::aes(x=log10GenesPerUMI, color=sample_id, fill=sample_id)) +
    ggplot2::geom_density(alpha=0.2) + ggplot2::theme_classic() +
    ggplot2::geom_vline(xintercept=0.85, linetype="dotted") +
    ggplot2::labs(title="log10(Genes per UMI)  (dotted = 0.85 suggestion)")
  ggplot2::ggsave(paste0(pfx_rm,"_log10GenesPerUMI_QC.pdf"), width=10, height=8)

  fs1 <- Seurat::FeatureScatter(data_seurat,"nCount_RNA","percent.mt",  group.by="sample_id") +
    ggplot2::geom_hline(yintercept=c(S$mito_hi,S$mito_lo), linetype="dotted") +
    ggplot2::scale_colour_manual(values=.okabe_pal)
  fs2 <- Seurat::FeatureScatter(data_seurat,"nCount_RNA","nFeature_RNA", group.by="sample_id") +
    ggplot2::geom_hline(yintercept=c(S$nGene_lo,S$nGene_hi), linetype="dotted") +
    ggplot2::scale_colour_manual(values=.okabe_pal)
  fs3 <- Seurat::FeatureScatter(data_seurat,"percent.ribo","percent.mt", group.by="sample_id") +
    ggplot2::geom_vline(xintercept=S$ribo_lo, linetype="dotted") +
    ggplot2::geom_hline(yintercept=c(S$mito_lo,S$mito_hi), linetype="dotted") +
    ggplot2::scale_colour_manual(values=.okabe_pal)
  ggplot2::ggsave(paste0(pfx_rm,"_all_QC_final.pdf"),
                  patchwork::wrap_plots(fs1, fs2, fs3, nrow = 1), width=20, height=8)

  log_info("QC plots saved to: ", dirs$ribomito)

  # ---- 7. Save ------------------------------------------------------------
  utils::write.csv(data_seurat[[]], paste0(pfx_csv,"_metadata_annotated.csv"), row.names=TRUE)
  save(meta_original, counts, data_seurat, aggr_info,
       file = paste0(pfx_R,"_seurat_annotated.RData"))
  saveRDS(data_seurat, file.path(dirs$scratch, "seurat_annotated.rds"))

  log_info("Step 02a complete. Output: ", qc_root)
  message("")
  message("  Inspect QC plots in: ", dirs$ribomito)
  message("  Then run ciri_step02_filter() with your real thresholds.")
  invisible(qc_root)
}


# =============================================================================
# ciri_step02_filter
# =============================================================================

#' Filter cells and build a Monocle3 CDS
#'
#' Applies QC thresholds to the annotated Seurat object produced by
#' \code{\link{ciri_step02_annotate}}, saves post-filter QC plots, then
#' builds a Monocle3 \code{cell_data_set} at 7 clustering resolutions so
#' you can inspect the knee plot and choose the best.
#'
#' @section Three-run workflow:
#' \enumerate{
#'   \item \strong{Run 1} — set filter thresholds, use default \code{num_dim = 30}.
#'         Inspect the knee plot (\code{_variance_knee_plot_dim30.png}).
#'   \item \strong{Run 2} — if the knee plot shows the curve has not flattened,
#'         re-run with a higher \code{num_dim}. Compare the 7 UMAP PDFs and pick
#'         a resolution.
#'   \item \strong{Run 3} — set \code{chosen_cds} to your preferred resolution.
#'         The final CDS is saved to \code{to_scratch/cds_final.rds}.
#' }
#'
#' @param scratch_dir  Folder containing \code{seurat_annotated.rds} from
#'                     \code{ciri_step02_annotate()}. Default \code{"scratch"}.
#' @param output_root  Top-level output directory. Default \code{"Output"}.
#' @param sample       Experiment name.
#' @param mito_lo      Lower mito bound: \code{percent.mt >= mito_lo}. Default \code{0}.
#' @param mito_hi      Upper mito bound: \code{percent.mt < mito_hi}. Default \code{10}.
#' @param ribo_lo      Lower ribo bound: \code{percent.ribo >= ribo_lo}. Default \code{1}.
#' @param nGene_lo     Min genes: \code{nFeature_RNA > nGene_lo}. Default \code{300}.
#' @param nGene_hi     Max genes: \code{nFeature_RNA < nGene_hi}. Default \code{7000}.
#' @param nUMI_lo      Min UMIs: \code{nCount_RNA > nUMI_lo}. Default \code{100}.
#' @param num_dim      PCA dimensions for Monocle3 \code{preprocess_cds()}.
#'                     Default \code{30}. Increase if the knee plot shows the
#'                     variance curve has not flattened.
#' @param chosen_cds   Which of the 7 clusterings to finalise as
#'                     \code{cds_final.rds}. One of \code{"cds_1"} through
#'                     \code{"cds_7"} (see Details). Default \code{"cds_3"}.
#'                     All 7 are always computed and saved regardless.
#'
#' @details
#' The 7 clustering resolutions, ordered from coarsest to finest:
#' \tabular{lll}{
#'   \strong{Name} \tab \strong{Resolution} \tab \strong{k} \cr
#'   cds_6 \tab 1e-5   \tab default \cr
#'   cds_7 \tab 2.5e-5 \tab default \cr
#'   cds_2 \tab 1e-4   \tab default \cr
#'   cds_4 \tab 2e-4   \tab default \cr
#'   cds_3 \tab 2e-4   \tab 15      \cr
#'   cds_5 \tab 5e-4   \tab default \cr
#'   cds_1 \tab 1e-2   \tab default \cr
#' }
#'
#' @return Invisibly returns the path to the Monocle3 output folder.
#'
#' @section Output structure:
#' \preformatted{
#' Output/QC/<date>/
#'   csv/        <date>_post_QCcounts.csv
#'               <date>_cell_metadata_filtered_<sample>.csv
#'   filtering/  <date>_all_QC_after_FILTER.pdf
#'   R_objects/  <date>_seurat_data_filtered.RData
#' Output/monocle/<date>/
#'   umap/       <date>_variance_knee_plot_dim<N>.png
#'               <date>_cds_1 to cds_7 UMAP PDFs
#'               <date>_UMAP.pdf  (chosen resolution)
#'   R_objects/  <date>_monocle_cds_clustered.RData
#'               <date>_monocle_cds_clustered_FINAL.RData
#'   to_scratch/ cds_final.rds
#' }
#'
#' @examples
#' \dontrun{
#' # Run 1 — apply thresholds, check knee plot
#' ciri_step02_filter(sample = "AB011", mito_hi = 10, nGene_lo = 300)
#'
#' # Run 2 — more PCA dims if needed
#' ciri_step02_filter(sample = "AB011", num_dim = 50, mito_hi = 10, nGene_lo = 300)
#'
#' # Run 3 — finalise chosen resolution
#' ciri_step02_filter(sample = "AB011", chosen_cds = "cds_3",
#'                    num_dim = 30, mito_hi = 10, nGene_lo = 300)
#' }
#' @export
ciri_step02_filter <- function(scratch_dir  = "scratch",
                                output_root  = "Output",
                                sample       = "CIRI",
                                mito_lo      = 0,
                                mito_hi      = 10,
                                ribo_lo      = 1,
                                nGene_lo     = 300,
                                nGene_hi     = 7000,
                                nUMI_lo      = 100,
                                num_dim      = 30L,
                                chosen_cds   = "cds_3") {

  .check_cran_pkgs("Seurat", "SeuratObject")
  .check_bioc_pkgs("monocle3")

  rds_path <- file.path(scratch_dir, "seurat_annotated.rds")
  assert_file(rds_path,
              hint = "Run ciri_step02_annotate() and promote to scratch/ first.")

  if (!chosen_cds %in% paste0("cds_", 1:7))
    stop("chosen_cds must be one of cds_1 … cds_7", call. = FALSE)

  # Output dirs
  d        <- format(Sys.Date(), "%y%m%d")
  qc_root  <- file.path(output_root, "QC",      d)
  mon_root <- file.path(output_root, "monocle",  d)
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

  F <- list(mito_lo=mito_lo, mito_hi=mito_hi, ribo_lo=ribo_lo,
            nGene_lo=nGene_lo, nGene_hi=nGene_hi, nUMI_lo=nUMI_lo)

  step_banner("02b", "Filter + Monocle3",
    inputs  = rds_path,
    outputs = c(qc_root, mon_root)
  )
  message(sprintf("  Thresholds: mito [%s,%s)  ribo>=%s  nGene (%s,%s)  nUMI>%s",
                  F$mito_lo, F$mito_hi, F$ribo_lo, F$nGene_lo, F$nGene_hi, F$nUMI_lo))
  message(sprintf("  num_dim: %d  |  chosen_cds: %s", num_dim, chosen_cds))

  # ---- 1. Load Seurat -----------------------------------------------------
  log_info("Loading seurat_annotated.rds ...")
  data_seurat  <- readRDS(rds_path)
  log_info("Loaded: ", ncol(data_seurat), " cells")
  pre_QCcounts <- data_seurat[[]] |>
    dplyr::group_by(sample_id) |>
    dplyr::summarise(pre_QC_count = dplyr::n(), .groups = "drop")

  # ---- 2. Filter ----------------------------------------------------------
  log_info("Applying filters ...")
  filtered_seurat <- base::subset(
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

  post_QCcounts <- metadata |>
    dplyr::group_by(sample_id) |>
    dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
    dplyr::mutate(
      threshold_nGene = paste0(F$nGene_lo, " - ", F$nGene_hi),
      threshold_nUMI  = paste0("> ",  F$nUMI_lo),
      threshold_mito  = paste0("[",   F$mito_lo, ", ", F$mito_hi, ")%"),
      threshold_ribo  = paste0(">= ", F$ribo_lo, "%")
    ) |>
    dplyr::left_join(pre_QCcounts, by = "sample_id") |>
    dplyr::mutate(cells_lost = pre_QC_count - count)

  utils::write.csv(post_QCcounts,
                   paste0(pfx_qc_csv,"_post_QCcounts.csv"), row.names = TRUE)
  utils::write.csv(metadata,
                   paste0(pfx_qc_csv,"_cell_metadata_filtered_",sample,".csv"),
                   row.names = TRUE)
  print(post_QCcounts)

  fs1 <- Seurat::FeatureScatter(filtered_seurat,"nCount_RNA","percent.mt",  group.by="sample_id") +
    ggplot2::geom_vline(xintercept=F$nUMI_lo, linetype="dotted") +
    ggplot2::geom_hline(yintercept=c(F$mito_lo,F$mito_hi), linetype="dotted") +
    ggplot2::scale_colour_manual(values=.pal15)
  fs2 <- Seurat::FeatureScatter(filtered_seurat,"nCount_RNA","nFeature_RNA", group.by="sample_id") +
    ggplot2::geom_vline(xintercept=F$nUMI_lo, linetype="dotted") +
    ggplot2::geom_hline(yintercept=c(F$nGene_lo,F$nGene_hi), linetype="dotted") +
    ggplot2::scale_colour_manual(values=.pal15)
  fs3 <- Seurat::FeatureScatter(filtered_seurat,"percent.ribo","percent.mt", group.by="sample_id") +
    ggplot2::geom_vline(xintercept=4, linetype="dotted") +
    ggplot2::geom_hline(yintercept=c(F$mito_lo,F$mito_hi), linetype="dotted") +
    ggplot2::scale_colour_manual(values=.pal15)
  ggplot2::ggsave(paste0(pfx_filter,"_all_QC_after_FILTER.pdf"),
                  patchwork::wrap_plots(fs1, fs2, fs3, nrow = 1), width=20, height=8)

  save(metadata, counts, filtered_seurat,
       file = paste0(pfx_qc_R,"_seurat_data_filtered.RData"))

  # ---- 3. Monocle3 CDS ----------------------------------------------------
  log_info("Building Monocle3 CDS ...")
  gene_annotation <- data.frame(gene_short_name = rownames(counts),
                                 row.names       = rownames(counts))
  cds <- monocle3::new_cell_data_set(counts,
                                      cell_metadata = metadata,
                                      gene_metadata = gene_annotation)

  log_info(sprintf("preprocess_cds(num_dim = %d) ...", num_dim))
  cds <- monocle3::preprocess_cds(cds, num_dim = num_dim)

  knee_file <- paste0(pfx_umap,"_variance_knee_plot_dim", num_dim, ".png")
  grDevices::png(knee_file, width=1200, height=800, res=100)
  print(monocle3::plot_pc_variance_explained(cds))
  grDevices::dev.off()

  message("")
  message("  ★ INSPECT KNEE PLOT: ", knee_file)
  message("    Curve should flatten before dim ", num_dim, ".")
  message("    If not, re-run with a larger num_dim, e.g. num_dim = 50")
  message("")

  set.seed(678686)
  cds <- monocle3::reduce_dimension(cds, reduction_method = "tSNE")
  cds <- monocle3::reduce_dimension(cds, umap.min_dist = 0.2,
                                     reduction_method = "UMAP")

  # ---- 4. Cluster at 7 resolutions ----------------------------------------
  log_info("Clustering at 7 resolutions ...")
  res_cfg <- list(
    cds_1 = list(res = 1e-2,   k = NULL, label = "resolution_1e_2"),
    cds_2 = list(res = 1e-4,   k = NULL, label = "resolution_1e_4"),
    cds_3 = list(res = 2e-4,   k = 15,   label = "resolution_2e_4_K15"),
    cds_4 = list(res = 2e-4,   k = NULL, label = "resolution_2e_4"),
    cds_5 = list(res = 5e-4,   k = NULL, label = "resolution_5e_4"),
    cds_6 = list(res = 1e-5,   k = NULL, label = "resolution_1e_5"),
    cds_7 = list(res = 2.5e-5, k = NULL, label = "resolution_25e_6")
  )

  cds_list <- list()
  for (nm in names(res_cfg)) {
    cfg  <- res_cfg[[nm]]
    args <- list(cds, resolution = cfg$res)
    if (!is.null(cfg$k)) args$k <- cfg$k
    clust <- do.call(monocle3::cluster_cells, args)
    SummarizedExperiment::colData(clust)$monocle_clusters <-
      as.character(monocle3::clusters(clust))
    cds_list[[nm]] <- clust

    p_umap <- monocle3::plot_cells(clust, color_cells_by="cluster",
                                    group_cells_by="cluster",
                                    group_label_size=7, cell_size=0.8)
    ggplot2::ggsave(paste0(pfx_umap,"_",nm,"_monocle_clusters_",cfg$label,".pdf"),
                    p_umap + ggplot2::facet_grid(~sample_id), width=20, height=4)
    ggplot2::ggsave(paste0(pfx_umap,"_",nm,"_monocle_clusters_",cfg$label,"_single.pdf"),
                    p_umap, width=4, height=4)
    log_info("  ", nm, " (res=", cfg$res,
             if (!is.null(cfg$k)) paste0(", k=",cfg$k) else "",
             ") — ", length(unique(monocle3::clusters(clust))), " clusters")
  }

  list2env(cds_list, envir = environment())   # cds_1 … cds_7 in local scope
  save(metadata, counts, cds, cds_1, cds_2, cds_3, cds_4, cds_5, cds_6, cds_7,
       file = paste0(pfx_mon_R,"_monocle_cds_clustered.RData"))
  saveRDS(cds_list, file.path(dirs$mon_scr,"cds_all_resolutions.rds"))

  # ---- 5. Finalise chosen CDS ---------------------------------------------
  log_info("Finalising chosen CDS: ", chosen_cds)
  cds_final <- monocle3::learn_graph(cds_list[[chosen_cds]])
  cds       <- cds_final

  pbuild <- monocle3::plot_cells(cds, color_cells_by="monocle_clusters",
                                  label_cell_groups=TRUE, show_trajectory_graph=FALSE,
                                  label_leaves=FALSE, label_branch_points=FALSE,
                                  graph_label_size=3, group_label_size=4, cell_size=1)
  SummarizedExperiment::colData(cds)$UMAP_1_monocle <- pbuild$data$data_dim_1
  SummarizedExperiment::colData(cds)$UMAP_2_monocle <- pbuild$data$data_dim_2

  p_final <- monocle3::plot_cells(cds, color_cells_by="cluster",
                                   group_cells_by="cluster", group_label_size=7,
                                   cell_size=0.8, show_trajectory_graph=FALSE)
  ggplot2::ggsave(paste0(pfx_umap,"_UMAP.pdf"),
                  p_final + ggplot2::facet_grid(~sample_id), width=10, height=5)
  ggplot2::ggsave(paste0(pfx_umap,"_UMAP_single.pdf"), p_final, width=5, height=5)

  cell_metadata <- as.data.frame(cds@colData)
  save(cell_metadata, counts, cds,
       file = paste0(pfx_mon_R,"_monocle_cds_clustered_FINAL.RData"))
  saveRDS(cds, file.path(dirs$mon_scr, "cds_final.rds"))

  log_info("Step 02b complete.")
  log_info("QC output      : ", qc_root)
  log_info("Monocle output : ", mon_root)
  message("")
  message("  Compare 7 UMAP resolutions in: ", dirs$mon_umap)
  message("  _cds_6_ res=1e-5  (fewest) … _cds_1_ res=1e-2 (most)")
  message("  Then re-run with chosen_cds = 'cds_N' to save final CDS.")
  message("  Final CDS: ", file.path(dirs$mon_scr, "cds_final.rds"))

  invisible(mon_root)
}
