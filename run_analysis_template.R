# =============================================================================
# CIRI Pipeline — Analysis Template
# =============================================================================
# Copy this file and rename it for each experiment, e.g.:
#   cp run_analysis_template.R run_analysis_AB011.R
#
# Run each block from the RStudio console one at a time.
# Inspect the output after each step before proceeding.
#
# Working directory should be the CIRI_pipeline/ repo root.
# (RStudio: Session → Set Working Directory → To Project Directory)
# =============================================================================

# ---- Configuration ----------------------------------------------------------
DATA_DIR    <- "/path/to/your/data"           # folder with .h5 and guides.csv
MATRIX      <- "filtered_feature_bc_matrix.h5"
SAMPLE      <- "AB011"                        # used in all output folder names
OUTPUT_ROOT <- "Output"                       # relative to project root
STEPS_DIR   <- "steps"

# Helper: build the Rscript call
run_step <- function(script, ...) {
  args <- c(file.path(STEPS_DIR, script), ...)
  message("\n>>> Running: ", paste(args, collapse = " "), "\n")
  system2("Rscript", args)
}

# =============================================================================
# STEP 00 — Download Ensembl gene reference (run once per project)
# =============================================================================
run_step("step_00_download_gene_ref.R",
  "--output_root", OUTPUT_ROOT,
  "--sample",      SAMPLE
)
# ↑ Inspect: Output/<date>_step00_geneRef_AB011/csv/ensembl_protein_coding_genes.csv

# =============================================================================
# STEP 01 — Perturbation Assignment
# =============================================================================
run_step("step_01_perturbation_assignment.R",
  "--data_dir",    DATA_DIR,
  "--matrix",      MATRIX,
  "--sample",      SAMPLE,
  "--output_root", OUTPUT_ROOT,
  "--strategy",    "1",       # 1 = single guide, 2 = dual guides
  "--threshold_a", "-1",      # -1 = auto-detect; set manually if KDE looks wrong
  "--threshold_i", "-1"
)
# ↑ Inspect before continuing:
#   Output/<date>_step01_assignment_AB011/plots/
#   - total_umixguide.pdf
#   - fixed_vs_variable_scatter.pdf
#   - single_guide_plots/
#   Output/<date>_step01_assignment_AB011/stats/assignment_summary.txt

# =============================================================================
# STEP 02 — Annotation & Filtering
# =============================================================================
run_step("step_02_anno_filter.R",
  "--data_dir",    DATA_DIR,
  "--matrix",      MATRIX,
  "--sample",      SAMPLE,
  "--output_root", OUTPUT_ROOT,
  "--min_genes",   "250",     # raise if you want stricter cell QC
  "--remove_mt",   "TRUE",    # removes ^MT- genes (protein-coding MT fix)
  "--remove_rb",   "TRUE"     # removes ^RPS|^RPL genes
)
# ↑ Inspect before continuing:
#   Output/<date>_step02_filter_AB011/plots/
#   - RiboMito_pre_filter.pdf
#   - RiboMito_post_cell_filter.pdf
#   Output/<date>_step02_filter_AB011/stats/filter_summary.txt

# =============================================================================
# STEP 03 — Load & Preprocess (Monocle3, UMAP, Clustering)
# =============================================================================
run_step("step_03_load.R",
  "--sample",      SAMPLE,
  "--output_root", OUTPUT_ROOT,
  "--resolution",  "5e-5"     # tune after looking at umap_clusters.pdf
                               # more clusters → increase; fewer → decrease
)
# ↑ Inspect before continuing:
#   Output/<date>_step03_load_AB011/plots/umap_clusters.pdf
#   Output/<date>_step03_load_AB011/stats/cluster_summary.txt
#   → Decide which cluster(s) are your target lineage for Steps 05–06

# =============================================================================
# STEP 04 — Target Validation (knockdown / activation efficiency)
# =============================================================================
run_step("step_04_target_validation.R",
  "--sample",      SAMPLE,
  "--output_root", OUTPUT_ROOT,
  "--control_a",   "NTCa-NA",  # adjust to match your non-targeting control name
  "--control_i",   "NTCa-NA"
)
# ↑ Inspect:
#   Output/<date>_step04_validation_AB011/plots/validation_violin_<gene>.pdf
#   Output/<date>_step04_validation_AB011/csv/validation_summary.csv

# =============================================================================
# STEP 05 — Cluster Enrichment Analysis
# =============================================================================
TARGET_CLUSTER <- "5"   # ← from umap_clusters.pdf (Step 03); use "3-4" for multiple

run_step("step_05_cluster_enrichment.R",
  "--sample",      SAMPLE,
  "--output_root", OUTPUT_ROOT,
  "--clusters",    TARGET_CLUSTER,
  "--control",     "NTCa-NA",
  "--min_cells",   "10"
)
# ↑ Inspect:
#   Output/<date>_step05_enrichment_AB011/plots/heatmap_enrichment_Group_5.pdf
#   Output/<date>_step05_enrichment_AB011/plots/scatter_enrichment_Group_5.pdf

# =============================================================================
# STEP 06 — Subclustering & Trajectory
# =============================================================================
GROUP     <- "muscle"   # ← short name for this lineage
ROOT_GENE <- "SOX2"     # ← marker gene highest at the start of differentiation

run_step("step_06_subclusters_trajectory.R",
  "--sample",      SAMPLE,
  "--output_root", OUTPUT_ROOT,
  "--clusters",    TARGET_CLUSTER,
  "--root_gene",   ROOT_GENE,
  "--group",       GROUP,
  "--resolution",  "1e-3"  # tune after looking at UMAP_muscle_subclusters.pdf
)
# ↑ Inspect before continuing:
#   Output/<date>_step06_trajectory_AB011/plots/
#   - UMAP_muscle_mainclusters.pdf
#   - UMAP_muscle_subclusters.pdf    ← tune --resolution if needed
#   - UMAP_muscle_pseudotime.pdf     ← check root placement looks correct

# =============================================================================
# STEP 07 — Pseudotime Statistics (KS Test)
# =============================================================================
run_step("step_07_pseudotime_stats.R",
  "--sample",         SAMPLE,
  "--output_root",    OUTPUT_ROOT,
  "--group",          GROUP,
  "--control",        "NTCa-NA",
  "--min_cells",      "8",
  "--ecdf_top_n",     "10",
  "--run_per_sample", "FALSE"
)
# ↑ Inspect:
#   Output/<date>_step07_pseudotime_AB011/plots/volcano_muscle.pdf
#   Output/<date>_step07_pseudotime_AB011/csv/ks_results_muscle.csv

# =============================================================================
# STEP 08 — Signature Scoring
# =============================================================================
run_step("step_08_signatures.R",
  "--sample",      SAMPLE,
  "--output_root", OUTPUT_ROOT,
  "--group",       GROUP
)
# ↑ Inspect:
#   Output/<date>_step08_signatures_AB011/plots/umap_<sig>.pdf
#   Output/<date>_step08_signatures_AB011/csv/signature_summary.csv
#   Output/<date>_step08_signatures_AB011/stats/signature_coverage.txt
#   → To add gene sets: edit SIGNATURES in steps/step_08_signatures.R
