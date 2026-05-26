# =============================================================================
# CIRI — Steps 02–08: exported pipeline functions
# =============================================================================

# ---- Step 02 ---------------------------------------------------------------

#' Filter and annotate the expression matrix
#'
#' Reads \code{annotation_data.csv} from \code{scratch/}, applies QC and gene
#' filters, and saves the annotated matrix to \code{to_scratch/}.
#'
#' @param data_dir    Directory containing the H5 matrix.
#' @param matrix      H5 matrix filename.
#' @param sample      Experiment name.
#' @param output_root Top-level output directory. Default \code{"Output"}.
#' @param scratch_dir Project-level scratch directory. Default \code{"scratch"}.
#' @param gene_ref    Path to \code{ensembl_protein_coding_genes.csv}.
#'                    If \code{NULL}, resolved from scratch/.
#' @param min_genes   Minimum detected genes per cell. Default \code{250}.
#' @param min_umis    Minimum total UMIs per gene. Default \code{3}.
#' @param remove_mt   Remove mitochondrial genes (\code{^MT-}). Default \code{TRUE}.
#' @param remove_rb   Remove ribosomal genes (\code{^RPS|^RPL}). Default \code{TRUE}.
#'
#' @return Invisibly returns the output folder path.
#' @examples
#' \dontrun{
#' ciri_step02_filter(
#'   data_dir = "/data/AB011",
#'   matrix   = "filtered_feature_bc_matrix.h5",
#'   sample   = "AB011"
#' )
#' }
#' @export
ciri_step02_filter <- function(data_dir,
                               matrix,
                               sample       = "CIRI",
                               output_root  = "Output",
                               scratch_dir  = "scratch",
                               gene_ref     = NULL,
                               min_genes    = 250L,
                               min_umis     = 3L,
                               remove_mt    = TRUE,
                               remove_rb    = TRUE) {

  .check_cran_pkgs("hdf5r", "data.table", "Seurat")
  .check_bioc_pkgs("monocle3")

  out      <- make_out_dirs(output_root, "step02_filter", sample)
  data_dir <- normalizePath(data_dir, mustWork = TRUE)

  ann_path <- file.path(scratch_dir, "annotation_data.csv")
  ref_path <- if (!is.null(gene_ref)) gene_ref else
    file.path(scratch_dir, "ensembl_protein_coding_genes.csv")

  step_banner("02", "Annotation & Filtering",
    inputs  = c(file.path(data_dir, matrix), ann_path, ref_path),
    outputs = c(out$plots, out$stats,
                file.path(out$to_scratch, "annotated_matrix.csv"))
  )

  assert_file(file.path(data_dir, matrix))
  assert_file(ann_path,  hint = "Run ciri_step01_assignment() then ciri_promote_scratch().")
  assert_file(ref_path,  hint = "Run ciri_step00_download_ref() then ciri_promote_scratch(), or pass gene_ref=.")

  log_info("Loading 10X data ...")
  data <- Seurat::Read10X_h5(file.path(data_dir, matrix), unique.features = TRUE)
  seu  <- Seurat::CreateSeuratObject(counts = data$`Gene Expression`, assay = "RNA")
  log_info("Loaded: ", ncol(seu), " cells, ", nrow(seu), " features")

  .mito_ribo_plot(seu, "Mito vs Ribo — all cells",
                  file.path(out$plots, "RiboMito_pre_filter.pdf"))

  n_before <- ncol(seu)
  seu <- Seurat::subset(seu, subset = nFeature_RNA >= min_genes)
  log_info("Cell filter (>= ", min_genes, "): ", n_before, " → ", ncol(seu))
  .mito_ribo_plot(seu, paste0("Mito vs Ribo — >= ", min_genes, " genes"),
                  file.path(out$plots, "RiboMito_post_cell_filter.pdf"))

  ref   <- read.csv(ref_path)
  if (!"mix" %in% names(ref))
    ref$mix <- ifelse(is.na(ref$hgnc_symbol) | ref$hgnc_symbol == "",
                      ref$ensembl_gene_id, ref$hgnc_symbol)
  mat   <- Seurat::GetAssayData(seu, assay = "RNA", layer = "counts")
  mat   <- mat[intersect(rownames(mat), ref$mix), ]
  log_info("Protein-coding filter: ", nrow(mat), " genes")

  n_mt <- 0L
  if (remove_mt) {
    mt <- grepl("^MT-", rownames(mat)); n_mt <- sum(mt)
    mat <- mat[!mt, ]; log_info("Removed ", n_mt, " MT genes → ", nrow(mat))
  }
  n_rb <- 0L
  if (remove_rb) {
    rb <- grepl("^RPS|^RPL", rownames(mat)); n_rb <- sum(rb)
    mat <- mat[!rb, ]; log_info("Removed ", n_rb, " RB genes → ", nrow(mat))
  }

  gene_umi <- Matrix::rowSums(mat)
  n_low    <- sum(gene_umi < min_umis)
  mat      <- mat[gene_umi >= min_umis, ]
  mat      <- mat[, Matrix::colSums(mat) > 0]
  log_info("Final matrix: ", nrow(mat), " genes × ", ncol(mat), " cells")

  assigned <- read.csv(ann_path)
  mat_df   <- as.data.frame(as.matrix(mat))
  colnames(mat_df) <- sub("\\.", "-", colnames(mat_df))
  shared   <- intersect(colnames(mat_df), assigned$cell_barcode)
  mat_df   <- mat_df[, shared]
  assigned <- dplyr::mutate(assigned, new_name = paste(cell_barcode, feature_a, feature_i, sep = "-"))
  name_map         <- stats::setNames(assigned$new_name, assigned$cell_barcode)
  colnames(mat_df) <- name_map[colnames(mat_df)]

  stats_txt <- c(paste("Input cells:", n_before), paste("After min_genes:", ncol(seu)),
                 paste("MT genes removed:", n_mt), paste("RB genes removed:", n_rb),
                 paste("Low-UMI genes removed:", n_low),
                 paste("Final genes:", nrow(mat_df)), paste("Final cells:", ncol(mat_df)))
  writeLines(stats_txt, file.path(out$stats, "filter_summary.txt"))
  message(paste(stats_txt, collapse = "\n"))

  data.table::fwrite(mat_df, file = file.path(out$to_scratch, "annotated_matrix.csv"),
                     row.names = TRUE)
  log_info("Step 02 complete. Output: ", out$out_dir)
  invisible(out$out_dir)
}

.mito_ribo_plot <- function(seu, title, outfile) {
  seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^MT-")
  seu[["percent.rb"]] <- Seurat::PercentageFeatureSet(seu, pattern = "^RPS|^RPL")
  lims <- range(seu@meta.data$nCount_RNA, na.rm = TRUE)
  p <- ggplot2::ggplot(seu@meta.data,
                       ggplot2::aes(x = percent.mt, y = percent.rb, colour = nCount_RNA)) +
    ggplot2::geom_point(alpha = 0.6, size = 1) + ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(x = "% Mito", y = "% Ribo", title = title) +
    ggplot2::scale_color_viridis_c(option = "turbo", trans = "log10",
                                   limits = lims, oob = scales::squish)
  ggplot2::ggsave(outfile, p, width = 10, height = 10)
}

# ---- Step 03 ---------------------------------------------------------------

#' Load annotated matrix into Monocle3 and preprocess
#'
#' Reads \code{annotated_matrix.csv} from \code{scratch/}, builds a Monocle3
#' \code{cell_data_set}, normalises, runs PCA, UMAP and Leiden clustering.
#'
#' @param sample      Experiment name.
#' @param output_root Top-level output directory. Default \code{"Output"}.
#' @param scratch_dir Project-level scratch directory. Default \code{"scratch"}.
#' @param resolution  Leiden clustering resolution. Default \code{5e-5}.
#' @param n_dims      PCA dimensions. Default \code{100}.
#' @param seed        Random seed. Default \code{1234597698}.
#'
#' @return Invisibly returns the output folder path.
#' @examples
#' \dontrun{
#' ciri_step03_load(sample = "AB011", resolution = 5e-5)
#' }
#' @export
ciri_step03_load <- function(sample       = "CIRI",
                             output_root  = "Output",
                             scratch_dir  = "scratch",
                             resolution   = 5e-5,
                             n_dims       = 100L,
                             seed         = 1234597698L) {

  .check_bioc_pkgs("monocle3")

  out <- make_out_dirs(output_root, "step03_load", sample)
  set.seed(seed)
  mat_path <- file.path(scratch_dir, "annotated_matrix.csv")

  step_banner("03", "Load & Preprocess",
    inputs  = mat_path,
    outputs = c(out$plots, out$R_objects,
                file.path(out$to_scratch, "processed_cds.RData"))
  )
  assert_file(mat_path, hint = "Run ciri_step02_filter() then ciri_promote_scratch().")

  log_info("Reading annotated matrix ...")
  exp           <- read.csv(mat_path, header = TRUE, check.names = FALSE)
  rownames(exp) <- exp[, 1]; exp[, 1] <- NULL
  log_info("Matrix: ", nrow(exp), " genes × ", ncol(exp), " cells")

  cell_meta <- data.frame(nomi = colnames(exp), row.names = colnames(exp)) |>
    tidyr::separate(nomi, into = c("cellID", "sample_id", "guide_a", "guide_i"),
                    sep = "-", remove = FALSE, fill = "right") |>
    dplyr::mutate(
      guide_a   = dplyr::na_if(guide_a, "NA"),
      guide_i   = dplyr::na_if(guide_i, "NA"),
      gene_a    = sapply(strsplit(guide_a, "_"), `[`, 1),
      gene_i    = sapply(strsplit(guide_i, "_"), `[`, 1),
      gene_comb = paste(gene_a, gene_i, sep = "-")
    )
  gene_meta <- data.frame(gene_short_name = rownames(exp), row.names = rownames(exp))

  log_info("Creating Monocle3 CDS ...")
  cds <- monocle3::new_cell_data_set(as.matrix(exp),
                                     cell_metadata = cell_meta,
                                     gene_metadata = gene_meta)
  log_info("Preprocessing ...")
  cds <- monocle3::preprocess_cds(cds, num_dim = n_dims)
  log_info("UMAP ...")
  cds <- monocle3::reduce_dimension(cds)
  log_info("Clustering (resolution = ", resolution, ") ...")
  cds <- monocle3::cluster_cells(cds, resolution = resolution)
  monocle3::colData(cds)$clusters <- monocle3::clusters(cds)

  n_clust <- length(unique(monocle3::clusters(cds)))
  log_info("Found ", n_clust, " clusters")

  .save_umap <- function(color_by, fname, label = FALSE) {
    pp <- monocle3::plot_cells(cds, color_cells_by = color_by,
                               show_trajectory_graph = FALSE,
                               label_cell_groups = label) +
      ggplot2::ggtitle(paste("UMAP —", color_by))
    ggplot2::ggsave(file.path(out$plots, fname), pp, width = 10, height = 8)
  }
  .save_umap("cluster", "umap_clusters.pdf", label = TRUE)
  .save_umap("sample",  "umap_sample.pdf")
  if ("gene_a" %in% names(monocle3::colData(cds))) .save_umap("gene_a", "umap_guide_a.pdf")
  if ("gene_i" %in% names(monocle3::colData(cds))) .save_umap("gene_i", "umap_guide_i.pdf")

  clust_tbl <- as.data.frame(table(monocle3::clusters(cds)))
  names(clust_tbl) <- c("cluster", "n_cells")
  utils::write.csv(clust_tbl, file.path(out$stats, "cluster_sizes.csv"), row.names = FALSE)

  rdata <- file.path(out$R_objects, "processed_cds.RData")
  save(cds, file = rdata)
  file.copy(rdata, file.path(out$to_scratch, "processed_cds.RData"), overwrite = TRUE)

  log_info("Step 03 complete. Output: ", out$out_dir)
  invisible(out$out_dir)
}

# ---- Step 04 ---------------------------------------------------------------

#' Validate target knockdown / activation efficiency
#'
#' Compares expression of each targeted gene in perturbed vs control cells.
#'
#' @param sample      Experiment name.
#' @param output_root Top-level output directory.
#' @param scratch_dir Project-level scratch directory.
#' @param control_a   Control \code{gene_comb} for CRISPRa. Default \code{"NTCa-NA"}.
#' @param control_i   Control \code{gene_comb} for CRISPRi. Default \code{"NTCa-NA"}.
#' @param min_cells   Minimum cells per group to generate a plot. Default \code{5}.
#'
#' @return Invisibly returns the output folder path.
#' @examples
#' \dontrun{
#' ciri_step04_validation(sample = "AB011", control_a = "NTCa-NA")
#' }
#' @export
ciri_step04_validation <- function(sample       = "CIRI",
                                   output_root  = "Output",
                                   scratch_dir  = "scratch",
                                   control_a    = "NTCa-NA",
                                   control_i    = "NTCa-NA",
                                   min_cells    = 5L) {

  .check_bioc_pkgs("monocle3")

  out      <- make_out_dirs(output_root, "step04_validation", sample)
  cds_path <- file.path(scratch_dir, "processed_cds.RData")

  step_banner("04", "Target Validation", inputs = cds_path,
              outputs = c(out$plots, out$csv))
  assert_file(cds_path, hint = "Run ciri_step03_load() then ciri_promote_scratch().")

  cds      <- load_rdata(cds_path)
  norm_mat <- monocle3::normalized_counts(cds)

  cell_df <- data.frame(nomi = colnames(norm_mat), row.names = colnames(norm_mat)) |>
    tidyr::separate(nomi, into = c("cellID", "samp", "guide_a", "guide_i"),
                    sep = "-", remove = FALSE, fill = "right") |>
    dplyr::mutate(guide_a = dplyr::na_if(guide_a, "NA"),
                  guide_i = dplyr::na_if(guide_i, "NA"),
                  gene_a  = sapply(strsplit(guide_a, "_"), `[`, 1),
                  gene_i  = sapply(strsplit(guide_i, "_"), `[`, 1),
                  gene_comb = paste(gene_a, gene_i, sep = "-"))

  targets <- dplyr::bind_rows(
    data.frame(gene = unique(cell_df$gene_a[!is.na(cell_df$gene_a) & cell_df$gene_comb != control_a]),
               class = "CRISPRa"),
    data.frame(gene = unique(cell_df$gene_i[!is.na(cell_df$gene_i) & cell_df$gene_comb != control_i]),
               class = "CRISPRi")
  ) |> dplyr::filter(!is.na(gene)) |> dplyr::distinct()

  rows <- list()
  for (i in seq_len(nrow(targets))) {
    gene  <- targets$gene[i]; class <- targets$class[i]
    if (!gene %in% rownames(norm_mat)) { log_warn("Gene not in matrix: ", gene); next }
    expr      <- norm_mat[gene, ]
    perturbed <- if (class == "CRISPRa") cell_df$nomi[!is.na(cell_df$gene_a) & cell_df$gene_a == gene] \
                 else                    cell_df$nomi[!is.na(cell_df$gene_i) & cell_df$gene_i == gene]
    control   <- if (class == "CRISPRa") cell_df$nomi[cell_df$gene_comb == control_a] \
                 else                    cell_df$nomi[cell_df$gene_comb == control_i]
    if (length(perturbed) < min_cells || length(control) < min_cells) {
      log_warn("Too few cells for ", gene, " — skipped"); next
    }
    df <- dplyr::bind_rows(
      data.frame(expr = as.numeric(expr[perturbed]), group = paste0(gene, " perturbed")),
      data.frame(expr = as.numeric(expr[control]),   group = "Control"))
    ggplot2::ggsave(file.path(out$plots, paste0("validation_violin_", gene, ".pdf")),
      ggplot2::ggplot(df, ggplot2::aes(x = group, y = expr, fill = group)) +
        ggplot2::geom_violin(trim = FALSE, alpha = 0.7) +
        ggplot2::geom_boxplot(width = 0.1, fill = "white", alpha = 0.7, outlier.size = 0.5) +
        ggplot2::theme_minimal(base_size = 13) +
        ggplot2::labs(title = paste(class, "—", gene), x = NULL, y = "Norm. expression") +
        ggplot2::theme(legend.position = "none"),
      width = 5, height = 5)
    rows[[gene]] <- data.frame(gene = gene, class = class,
      n_perturbed = length(perturbed), n_control = length(control),
      median_perturbed = stats::median(as.numeric(expr[perturbed])),
      median_control   = stats::median(as.numeric(expr[control])),
      fold_change = stats::median(as.numeric(expr[perturbed])) /
                    (stats::median(as.numeric(expr[control])) + 1e-9))
  }
  if (length(rows)) {
    df_out <- dplyr::bind_rows(rows)
    utils::write.csv(df_out, file.path(out$csv, "validation_summary.csv"), row.names = FALSE)
    print(df_out)
  }
  log_info("Step 04 complete. Output: ", out$out_dir)
  invisible(out$out_dir)
}

# ---- Step 05 ---------------------------------------------------------------

#' Cluster enrichment analysis
#'
#' Calculates what fraction of each perturbation ends up in target cluster(s).
#'
#' @param clusters    Target cluster ID(s) as a character vector, e.g. \code{"5"}
#'                    or \code{c("3","4")}.
#' @param control     Control \code{gene_comb} label.
#' @param sample      Experiment name.
#' @param output_root Top-level output directory.
#' @param scratch_dir Project-level scratch directory.
#' @param min_cells   Minimum cells per perturbation to include. Default \code{10}.
#'
#' @return Invisibly returns the output folder path.
#' @examples
#' \dontrun{
#' ciri_step05_enrichment(clusters = "5", control = "NTCa-NA", sample = "AB011")
#' }
#' @export
ciri_step05_enrichment <- function(clusters,
                                   control,
                                   sample       = "CIRI",
                                   output_root  = "Output",
                                   scratch_dir  = "scratch",
                                   min_cells    = 10L) {

  .check_bioc_pkgs("monocle3")

  out         <- make_out_dirs(output_root, "step05_enrichment", sample)
  cds_path    <- file.path(scratch_dir, "processed_cds.RData")
  group_label <- paste("Group", paste(clusters, collapse = "_"), sep = "_")

  step_banner("05", "Cluster Enrichment", inputs = cds_path,
              outputs = c(out$csv, out$plots))
  assert_file(cds_path, hint = "Run ciri_step03_load() then ciri_promote_scratch().")

  cds <- load_rdata(cds_path)
  monocle3::colData(cds)$clusters <- monocle3::clusters(cds)
  df  <- as.data.frame(monocle3::colData(cds))

  .overview <- function(gene_col, fname) {
    if (!gene_col %in% names(df)) return(invisible(NULL))
    counts <- df |> dplyr::group_by(sample, clusters, .data[[gene_col]]) |>
      dplyr::summarise(count = dplyr::n(), .groups = "drop") |>
      dplyr::group_by(sample, clusters) |>
      dplyr::mutate(pct = count / sum(count) * 100) |> dplyr::ungroup()
    n_v <- length(unique(counts[[gene_col]]))
    p <- ggplot2::ggplot(counts, ggplot2::aes(x = factor(sample), y = pct,
                                               fill = .data[[gene_col]])) +
      ggplot2::geom_bar(stat = "identity", colour = "white", linewidth = 0.2) +
      ggplot2::scale_fill_manual(values = viridis::turbo(n_v)) +
      ggplot2::facet_grid(~clusters) + ggplot2::theme_minimal() +
      ggplot2::labs(title = paste(gene_col, "across clusters"), x = "Sample", y = "% cells")
    ggplot2::ggsave(file.path(out$plots, fname), p, width = 20, height = 10)
  }
  .overview("gene_a", "gene_a_in_clusters.pdf")
  .overview("gene_i", "gene_i_in_clusters.pdf")

  enrichment <- df |> dplyr::mutate(in_target = clusters %in% clusters) |>
    dplyr::group_by(sample, gene_comb) |>
    dplyr::summarise(n_total = dplyr::n(), n_target = sum(in_target),
                     pct_target = n_target / n_total * 100, .groups = "drop") |>
    dplyr::filter(n_total >= min_cells)

  utils::write.csv(enrichment,
    file.path(out$csv, paste0("cluster_enrichment_", group_label, ".csv")), row.names = FALSE)

  samples <- unique(enrichment$sample)
  if (length(samples) >= 2) {
    wide <- enrichment |> dplyr::select(sample, gene_comb, pct_target) |>
      tidyr::pivot_wider(names_from = sample, values_from = pct_target, values_fill = NA)
    complete <- wide |> dplyr::filter(rowSums(is.na(dplyr::across(-gene_comb))) == 0) |>
      dplyr::pull(gene_comb)
    heat_df <- dplyr::filter(enrichment, gene_comb %in% complete)
    ggplot2::ggsave(file.path(out$plots, paste0("heatmap_enrichment_", group_label, ".pdf")),
      ggplot2::ggplot(heat_df, ggplot2::aes(x = factor(sample), y = gene_comb, fill = pct_target)) +
        ggplot2::geom_tile(colour = "white") +
        ggplot2::scale_fill_viridis_c(option = "plasma", name = "% in target") +
        ggplot2::theme_minimal(base_size = 10) +
        ggplot2::labs(title = paste("Enrichment — cluster(s)", paste(clusters, collapse = "+")),
                      x = "Sample", y = "Perturbation"),
      width = 8, height = max(5, length(complete) * 0.25))
  }
  log_info("Step 05 complete. Output: ", out$out_dir)
  invisible(out$out_dir)
}

# ---- Step 06 ---------------------------------------------------------------

#' Subcluster and learn pseudotime trajectory
#'
#' @param clusters    Cluster(s) to subset, character vector.
#' @param root_gene   Gene marking the trajectory root.
#' @param group       Short lineage name used in filenames.
#' @param sample      Experiment name.
#' @param output_root Top-level output directory.
#' @param scratch_dir Project-level scratch directory.
#' @param resolution  Re-clustering resolution within the subset. Default \code{1e-3}.
#' @param n_dims      PCA dimensions for subset. Default \code{50}.
#' @param seed        Random seed. Default \code{42}.
#'
#' @return Invisibly returns the output folder path.
#' @examples
#' \dontrun{
#' ciri_step06_trajectory(
#'   clusters  = "5",
#'   root_gene = "SOX2",
#'   group     = "muscle",
#'   sample    = "AB011"
#' )
#' }
#' @export
ciri_step06_trajectory <- function(clusters,
                                   root_gene,
                                   group,
                                   sample       = "CIRI",
                                   output_root  = "Output",
                                   scratch_dir  = "scratch",
                                   resolution   = 1e-3,
                                   n_dims       = 50L,
                                   seed         = 42L) {

  .check_bioc_pkgs("monocle3")
  .check_cran_pkgs("igraph")

  out      <- make_out_dirs(output_root, "step06_trajectory", sample)
  cds_path <- file.path(scratch_dir, "processed_cds.RData")
  set.seed(seed)

  step_banner("06", "Subclustering & Trajectory", inputs = cds_path,
              outputs = c(out$plots, out$csv, out$R_objects,
                          file.path(out$to_scratch, paste0("processed_cds_", group, ".RData")),
                          file.path(out$to_scratch, paste0("pseudotime_", group, ".csv"))))
  assert_file(cds_path, hint = "Run ciri_step03_load() then ciri_promote_scratch().")

  cds <- load_rdata(cds_path)
  monocle3::colData(cds)$clusters <- monocle3::clusters(cds)
  cds_sub <- cds[, monocle3::colData(cds)$clusters %in% clusters]
  monocle3::colData(cds_sub)$clusters_main <- monocle3::colData(cds_sub)$clusters
  log_info("Subset: ", ncol(cds_sub), " cells")

  cds_sub <- monocle3::preprocess_cds(cds_sub, num_dim = n_dims)
  cds_sub <- monocle3::reduce_dimension(cds_sub, reduction_method = "UMAP")
  cds_sub <- monocle3::cluster_cells(cds_sub, resolution = resolution, random_seed = seed)
  monocle3::colData(cds_sub)$clusters_sub <- monocle3::clusters(cds_sub)

  .su <- function(color_by, fname, label = FALSE, traj = FALSE) {
    pp <- monocle3::plot_cells(cds_sub, color_cells_by = color_by,
                               show_trajectory_graph = traj,
                               label_cell_groups = label) +
      ggplot2::ggtitle(paste("UMAP —", color_by, "|", group))
    ggplot2::ggsave(file.path(out$plots, fname), pp, width = 6, height = 5)
  }
  .su("clusters_main", paste0("UMAP_", group, "_mainclusters.pdf"),  label = TRUE)
  .su("clusters_sub",  paste0("UMAP_", group, "_subclusters.pdf"),   label = TRUE)

  cds_sub <- monocle3::learn_graph(cds_sub)

  root_node <- NULL
  if (root_gene %in% rownames(cds_sub)) {
    root_expr   <- as.numeric(monocle3::normalized_counts(cds_sub)[root_gene, ])
    graph_nodes <- igraph::V(monocle3::principal_graph(cds_sub)[["UMAP"]])$name
    closest     <- cds_sub@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
    mean_node   <- sapply(graph_nodes, function(nd) {
      idx <- which(closest == nd); if (!length(idx)) return(-Inf); mean(root_expr[idx]) })
    root_node   <- graph_nodes[which.max(mean_node)]
    log_info("Root node (", root_gene, "): ", root_node)
  }
  cds_sub <- tryCatch(
    monocle3::order_cells(cds_sub, root_pr_nodes = root_node),
    error = function(e) { log_warn("order_cells failed, running without root."); monocle3::order_cells(cds_sub) })

  .su("pseudotime",   paste0("UMAP_", group, "_pseudotime.pdf"),  traj = TRUE)
  .su("clusters_sub", paste0("UMAP_", group, "_trajectory.pdf"),  traj = TRUE, label = TRUE)

  pt_df   <- data.frame(pseudotime = monocle3::pseudotime(cds_sub), row.names = colnames(cds_sub))
  pt_file <- paste0("pseudotime_", group, ".csv")
  utils::write.csv(pt_df, file.path(out$csv, pt_file))
  file.copy(file.path(out$csv, pt_file), file.path(out$to_scratch, pt_file), overwrite = TRUE)

  cds_file  <- paste0("processed_cds_", group, ".RData")
  rdata     <- file.path(out$R_objects, cds_file)
  save(cds_sub, file = rdata)
  file.copy(rdata, file.path(out$to_scratch, cds_file), overwrite = TRUE)

  log_info("Step 06 complete. Output: ", out$out_dir)
  invisible(out$out_dir)
}

# ---- Step 07 ---------------------------------------------------------------

#' Pseudotime KS test
#'
#' Compares pseudotime distributions of perturbations vs control using a
#' Kolmogorov-Smirnov test.
#'
#' @param group          Group name matching \code{\link{ciri_step06_trajectory}}.
#' @param control        Control \code{gene_comb} string.
#' @param sample         Experiment name.
#' @param output_root    Top-level output directory.
#' @param scratch_dir    Project-level scratch directory.
#' @param min_cells      Min cells per perturbation. Default \code{8}.
#' @param run_per_sample Also run per-sample. Default \code{FALSE}.
#' @param ecdf_top_n     Number of top hits to plot as ECDFs. Default \code{10}.
#'
#' @return Invisibly returns the output folder path.
#' @examples
#' \dontrun{
#' ciri_step07_pseudotime(group = "muscle", control = "NTCa-NA", sample = "AB011")
#' }
#' @export
ciri_step07_pseudotime <- function(group,
                                   control,
                                   sample         = "CIRI",
                                   output_root    = "Output",
                                   scratch_dir    = "scratch",
                                   min_cells      = 8L,
                                   run_per_sample = FALSE,
                                   ecdf_top_n     = 10L) {

  .check_bioc_pkgs("monocle3")

  out      <- make_out_dirs(output_root, "step07_pseudotime", sample)
  cds_path <- file.path(scratch_dir, paste0("processed_cds_", group, ".RData"))
  pt_path  <- file.path(scratch_dir, paste0("pseudotime_", group, ".csv"))

  step_banner("07", "Pseudotime Statistics (KS Test)",
              inputs = c(cds_path, pt_path), outputs = c(out$csv, out$plots))
  assert_file(cds_path, hint = "Run ciri_step06_trajectory() then ciri_promote_scratch().")
  assert_file(pt_path,  hint = "Run ciri_step06_trajectory() then ciri_promote_scratch().")

  cds     <- load_rdata(cds_path)
  pt      <- utils::read.csv(pt_path, row.names = 1); names(pt) <- "pseudotime"
  cell_df <- as.data.frame(monocle3::colData(cds))
  shared  <- intersect(rownames(pt), rownames(cell_df))
  cell_df <- cell_df[shared, ]; cell_df$pseudotime <- pt[shared, "pseudotime"]

  .run_ks <- function(df, ctrl, min_c) {
    ctrl_pt <- df$pseudotime[df$gene_comb == ctrl]
    ctrl_pt <- ctrl_pt[!is.na(ctrl_pt) & is.finite(ctrl_pt)]
    if (length(ctrl_pt) < min_c) stop("Control has < ", min_c, " finite pseudotime values.")
    rows <- lapply(unique(df$gene_comb), function(comb) {
      if (comb == ctrl) return(NULL)
      pert_pt <- df$pseudotime[df$gene_comb == comb]
      pert_pt <- pert_pt[!is.na(pert_pt) & is.finite(pert_pt)]
      if (length(pert_pt) < min_c) return(NULL)
      ks <- stats::ks.test(pert_pt, ctrl_pt)
      data.frame(gene_comb = comb, n_cells = length(pert_pt), ks_stat = ks$statistic,
                 p_value = ks$p.value,
                 direction = ifelse(stats::median(pert_pt) > stats::median(ctrl_pt), "faster", "slower"),
                 median_pert = stats::median(pert_pt), median_ctrl = stats::median(ctrl_pt))
    })
    res <- dplyr::bind_rows(rows)
    res$p_adj <- stats::p.adjust(res$p_value, method = "BH")
    dplyr::arrange(res, p_adj, dplyr::desc(ks_stat))
  }

  ks_res <- .run_ks(cell_df, control, min_cells)
  utils::write.csv(ks_res, file.path(out$csv, paste0("ks_results_", group, ".csv")), row.names = FALSE)
  n_sig <- sum(ks_res$p_adj < 0.05, na.rm = TRUE)
  log_info("Significant hits (p.adj<0.05): ", n_sig)

  ks_plot <- dplyr::mutate(ks_res, log10_padj = -log10(p_adj + 1e-300),
    label = ifelse(p_adj < 0.05 & dplyr::row_number() <= 20, gene_comb, NA_character_))
  ggplot2::ggsave(file.path(out$plots, paste0("volcano_", group, ".pdf")),
    ggplot2::ggplot(ks_plot, ggplot2::aes(x = ks_stat, y = log10_padj, colour = direction, label = label)) +
      ggplot2::geom_point(ggplot2::aes(size = n_cells), alpha = 0.7) +
      ggrepel::geom_text_repel(na.rm = TRUE, size = 3) +
      ggplot2::geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey50") +
      ggplot2::scale_colour_manual(values = c(faster = "#2166AC", slower = "#D6604D")) +
      ggplot2::theme_minimal(base_size = 13) +
      ggplot2::labs(title = paste("KS test —", group), x = "KS statistic",
                    y = "-log10(adj. p-value)", colour = "Direction", size = "n cells"),
    width = 10, height = 8)

  ctrl_pt <- cell_df$pseudotime[cell_df$gene_comb == control]
  ctrl_pt <- ctrl_pt[!is.na(ctrl_pt) & is.finite(ctrl_pt)]
  for (comb in utils::head(ks_res$gene_comb, ecdf_top_n)) {
    pert_pt <- cell_df$pseudotime[cell_df$gene_comb == comb]
    pert_pt <- pert_pt[!is.na(pert_pt) & is.finite(pert_pt)]
    ecdf_df <- dplyr::bind_rows(data.frame(pseudotime = pert_pt, group = comb),
                                 data.frame(pseudotime = ctrl_pt, group = control))
    safe    <- gsub("[^A-Za-z0-9_.-]", "_", comb)
    ggplot2::ggsave(file.path(out$plots, paste0("ecdf_", safe, "_", group, ".pdf")),
      ggplot2::ggplot(ecdf_df, ggplot2::aes(x = pseudotime, colour = group)) +
        ggplot2::stat_ecdf(linewidth = 1) + ggplot2::theme_minimal(base_size = 12) +
        ggplot2::scale_colour_manual(values = c("steelblue", "grey50"),
                                     breaks = c(comb, control)) +
        ggplot2::labs(title = paste("ECDF —", comb, "vs", control),
                      x = "Pseudotime", y = "Cumulative fraction", colour = NULL),
      width = 7, height = 5)
  }

  writeLines(c(paste("Group:", group), paste("Control:", control),
               paste("N tested:", nrow(ks_res)), paste("N significant:", n_sig),
               "Top 10:", utils::capture.output(print(utils::head(ks_res, 10)))),
             file.path(out$stats, paste0("ks_summary_", group, ".txt")))

  log_info("Step 07 complete. Output: ", out$out_dir)
  invisible(out$out_dir)
}

# ---- Step 08 ---------------------------------------------------------------

#' Score cells against gene set signatures
#'
#' @param group       Group name matching \code{\link{ciri_step06_trajectory}}.
#' @param sample      Experiment name.
#' @param output_root Top-level output directory.
#' @param scratch_dir Project-level scratch directory.
#' @param signatures  Named list of character vectors (gene symbols).
#'                    Defaults to built-in myogenic signatures if \code{NULL}.
#'
#' @return Invisibly returns the output folder path.
#' @examples
#' \dontrun{
#' ciri_step08_signatures(group = "muscle", sample = "AB011")
#' # custom signatures:
#' ciri_step08_signatures(
#'   group = "muscle", sample = "AB011",
#'   signatures = list(my_sig = c("MYOD1","MYOG","MYF5"))
#' )
#' }
#' @export
ciri_step08_signatures <- function(group,
                                   sample       = "CIRI",
                                   output_root  = "Output",
                                   scratch_dir  = "scratch",
                                   signatures   = NULL) {

  .check_bioc_pkgs("monocle3")

  out      <- make_out_dirs(output_root, "step08_signatures", sample)
  cds_path <- file.path(scratch_dir, paste0("processed_cds_", group, ".RData"))
  dir.create(file.path(out$csv, "signature_values"), showWarnings = FALSE)

  step_banner("08", "Signature Scoring", inputs = cds_path,
              outputs = c(out$plots, out$csv))
  assert_file(cds_path, hint = "Run ciri_step06_trajectory() then ciri_promote_scratch().")

  if (is.null(signatures)) signatures <- .default_signatures()

  cds     <- load_rdata(cds_path)
  present <- rownames(cds)
  rows    <- list()
  cov_log <- character(0)

  for (sig_name in names(signatures)) {
    genes   <- intersect(signatures[[sig_name]], present)
    n_tot   <- length(signatures[[sig_name]]); n_found <- length(genes)
    cov_log <- c(cov_log, sprintf("%-30s %d/%d (%.0f%%)", sig_name, n_found, n_tot, n_found/n_tot*100))
    log_info(utils::tail(cov_log, 1))
    if (n_found < 3) { log_warn("Skipping ", sig_name, " — < 3 genes."); next }

    norm_sub <- monocle3::normalized_counts(cds)[genes, , drop = FALSE]
    z_mat    <- t(scale(t(as.matrix(norm_sub))))
    scores   <- colMeans(z_mat, na.rm = TRUE)

    utils::write.csv(data.frame(score = scores, row.names = names(scores)),
                     file.path(out$csv, "signature_values", paste0(sig_name, ".csv")))
    monocle3::colData(cds)[[paste0("sig_", sig_name)]] <- scores[colnames(cds)]

    ggplot2::ggsave(file.path(out$plots, paste0("umap_", sig_name, ".pdf")),
      monocle3::plot_cells(cds, color_cells_by = paste0("sig_", sig_name),
                           show_trajectory_graph = FALSE, label_cell_groups = FALSE) +
        ggplot2::scale_color_viridis_c(option = "plasma", name = "Score") +
        ggplot2::ggtitle(paste("Signature:", sig_name)),
      width = 8, height = 6)

    if ("gene_comb" %in% names(monocle3::colData(cds))) {
      rows[[sig_name]] <- data.frame(gene_comb = monocle3::colData(cds)$gene_comb,
                                      score = scores[colnames(cds)]) |>
        dplyr::group_by(gene_comb) |>
        dplyr::summarise(mean_score = mean(score, na.rm = TRUE), n_cells = dplyr::n(),
                         .groups = "drop") |>
        dplyr::mutate(signature = sig_name)
    }
  }
  writeLines(cov_log, file.path(out$stats, "signature_coverage.txt"))
  if (length(rows)) utils::write.csv(dplyr::bind_rows(rows),
    file.path(out$csv, "signature_summary.csv"), row.names = FALSE)

  log_info("Step 08 complete. Output: ", out$out_dir)
  invisible(out$out_dir)
}

.default_signatures <- function() {
  list(
    cellcycle_s     = c("MCM5","PCNA","TYMS","FEN1","MCM2","MCM4","RRM1","UNG","GINS2",
                        "MCM6","CDCA7","DTL","PRIM1","UHRF1","HELLS","RFC2","RPA2","NASP",
                        "RAD51","CHEK1","ORC6","MCM3","MCM7","BRIP1","E2F8"),
    cellcycle_g2m   = c("HMGB2","CDK1","NUSAP1","UBE2C","BIRC5","TPX2","TOP2A","NDC80",
                        "CKS2","NUF2","CKS1B","MKI67","TMPO","CENPF","TACC3","SMC4",
                        "CCNB2","AURKB","BUB1","KIF11","CDC20","TTK"),
    sarcomere_core  = c("MYH3","MYH8","MYBPC1","MYBPC2","MYL1","MYL4","TNNI1","TNNI2",
                        "TNNT1","TNNT2","TNNC1","TNNC2","TPM1","TPM2","ACTA1","ACTC1","NEB","NEBL"),
    myogenic        = c("MYOD1","MYOG","MYF5","MYF6","PAX7","MEF2A","MEF2C","MEF2D",
                        "CDH15","NCAM1","DES","VIM","NES"),
    pluripotency    = c("POU5F1","SOX2","NANOG","KLF4","MYC","LIN28A","DPPA4","DPPA5",
                        "SALL4","UTF1","DNMT3L")
  )
}
