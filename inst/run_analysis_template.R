# =============================================================================
# CIRI Analysis Template
# =============================================================================
# Copy and rename this file per experiment, e.g.:
#   cp run_analysis_template.R run_analysis_AB011.R
#
# SETUP: either install the package from GitHub, or source the steps directly.
# See the two options below.
#
# Run each block from the RStudio console ONE AT A TIME.
# Inspect output after each step before promoting to scratch/.
# =============================================================================

# ---- Install options (choose one) ------------------------------------------

## Option A: install from GitHub (recommended, one-time)
# install.packages("devtools")
# devtools::install_github("yourusername/CIRI")
library(CIRI)

## Option B: source locally (no install needed, good for development)
# invisible(lapply(list.files("R", pattern="\\.R$", full.names=TRUE), source))

# ---- Configuration ----------------------------------------------------------
DATA_DIR    <- "/path/to/your/data"
MATRIX      <- "filtered_feature_bc_matrix.h5"
SAMPLE      <- "AB011"
OUTPUT_ROOT <- "Output"
SCRATCH_DIR <- "scratch"

# =============================================================================
# STEP 00 — Download Ensembl reference (run once per project)
# =============================================================================
ciri_step00_download_ref(
  output_root = OUTPUT_ROOT,
  sample      = SAMPLE
)
# Inspect: Output/<date>_step00_geneRef_AB011/csv/ensembl_protein_coding_genes.csv
# Then promote when happy:
ciri_promote_scratch("step00_geneRef", sample = SAMPLE,
                     output_root = OUTPUT_ROOT, scratch_dir = SCRATCH_DIR)

# =============================================================================
# STEP 01 — Perturbation Assignment
# =============================================================================
ciri_step01_assignment(
  data_dir    = DATA_DIR,
  matrix      = MATRIX,
  sample      = SAMPLE,
  output_root = OUTPUT_ROOT,
  strategy    = 1,       # 1 = single guide, 2 = dual guides
  threshold_a = -1,      # -1 = auto-detect via KDE valley
  threshold_i = -1
)
# ↑ Inspect before promoting:
#   Output/<date>_step01_assignment_AB011/plots/total_umixguide.pdf
#   Output/<date>_step01_assignment_AB011/plots/fixed_vs_variable_scatter.pdf
#   Output/<date>_step01_assignment_AB011/stats/assignment_summary.txt

ciri_promote_scratch("step01_assignment", sample = SAMPLE,
                     output_root = OUTPUT_ROOT, scratch_dir = SCRATCH_DIR)

# =============================================================================
# STEP 02 — Annotation & Filtering
# =============================================================================
ciri_step02_filter(
  data_dir    = DATA_DIR,
  matrix      = MATRIX,
  sample      = SAMPLE,
  output_root = OUTPUT_ROOT,
  scratch_dir = SCRATCH_DIR,
  min_genes   = 250,
  remove_mt   = TRUE,
  remove_rb   = TRUE
)
# ↑ Inspect before promoting:
#   Output/<date>_step02_filter_AB011/plots/RiboMito_pre_filter.pdf
#   Output/<date>_step02_filter_AB011/plots/RiboMito_post_cell_filter.pdf
#   Output/<date>_step02_filter_AB011/stats/filter_summary.txt

ciri_promote_scratch("step02_filter", sample = SAMPLE,
                     output_root = OUTPUT_ROOT, scratch_dir = SCRATCH_DIR)

# =============================================================================
# STEP 03 — Load & Preprocess (Monocle3, UMAP, Clustering)
# =============================================================================
ciri_step03_load(
  sample       = SAMPLE,
  output_root  = OUTPUT_ROOT,
  scratch_dir  = SCRATCH_DIR,
  resolution   = 5e-5     # tune after inspecting umap_clusters.pdf
)
# ↑ Inspect before promoting:
#   Output/<date>_step03_load_AB011/plots/umap_clusters.pdf
#   Output/<date>_step03_load_AB011/stats/cluster_sizes.csv
#   → decide which cluster is your target lineage (for steps 05-06)

ciri_promote_scratch("step03_load", sample = SAMPLE,
                     output_root = OUTPUT_ROOT, scratch_dir = SCRATCH_DIR)

# =============================================================================
# STEP 04 — Target Validation
# =============================================================================
ciri_step04_validation(
  sample      = SAMPLE,
  output_root = OUTPUT_ROOT,
  scratch_dir = SCRATCH_DIR,
  control_a   = "NTCa-NA",  # adjust to your NTC label
  control_i   = "NTCa-NA"
)
# ↑ Inspect:
#   Output/<date>_step04_validation_AB011/plots/validation_violin_<gene>.pdf
#   Output/<date>_step04_validation_AB011/csv/validation_summary.csv
# No promote needed — step 04 is a terminal analysis.

# =============================================================================
# STEP 05 — Cluster Enrichment
# =============================================================================
TARGET_CLUSTER <- c("5")  # from umap_clusters.pdf; e.g. c("3","4") for multiple

ciri_step05_enrichment(
  clusters    = TARGET_CLUSTER,
  control     = "NTCa-NA",
  sample      = SAMPLE,
  output_root = OUTPUT_ROOT,
  scratch_dir = SCRATCH_DIR,
  min_cells   = 10
)
# ↑ Inspect:
#   Output/<date>_step05_enrichment_AB011/plots/heatmap_enrichment_Group_5.pdf
# No promote needed — terminal analysis.

# =============================================================================
# STEP 06 — Subclustering & Trajectory
# =============================================================================
GROUP     <- "muscle"
ROOT_GENE <- "SOX2"

ciri_step06_trajectory(
  clusters    = TARGET_CLUSTER,
  root_gene   = ROOT_GENE,
  group       = GROUP,
  sample      = SAMPLE,
  output_root = OUTPUT_ROOT,
  scratch_dir = SCRATCH_DIR,
  resolution  = 1e-3    # tune after inspecting UMAP_muscle_subclusters.pdf
)
# ↑ Inspect before promoting:
#   Output/<date>_step06_trajectory_AB011/plots/UMAP_muscle_subclusters.pdf
#   Output/<date>_step06_trajectory_AB011/plots/UMAP_muscle_pseudotime.pdf
#   → check root placement looks correct

ciri_promote_scratch("step06_trajectory", sample = SAMPLE,
                     output_root = OUTPUT_ROOT, scratch_dir = SCRATCH_DIR)

# =============================================================================
# STEP 07 — Pseudotime Statistics (KS Test)
# =============================================================================
ciri_step07_pseudotime(
  group          = GROUP,
  control        = "NTCa-NA",
  sample         = SAMPLE,
  output_root    = OUTPUT_ROOT,
  scratch_dir    = SCRATCH_DIR,
  min_cells      = 8,
  ecdf_top_n     = 10,
  run_per_sample = FALSE
)
# ↑ Inspect:
#   Output/<date>_step07_pseudotime_AB011/plots/volcano_muscle.pdf
#   Output/<date>_step07_pseudotime_AB011/csv/ks_results_muscle.csv
# No promote needed — terminal analysis.

# =============================================================================
# STEP 08 — Signature Scoring
# =============================================================================
ciri_step08_signatures(
  group       = GROUP,
  sample      = SAMPLE,
  output_root = OUTPUT_ROOT,
  scratch_dir = SCRATCH_DIR
  # signatures = list(my_sig = c("MYOD1","MYOG")) # optional custom gene sets
)
# ↑ Inspect:
#   Output/<date>_step08_signatures_AB011/plots/umap_<sig>.pdf
#   Output/<date>_step08_signatures_AB011/csv/signature_summary.csv
#   Output/<date>_step08_signatures_AB011/stats/signature_coverage.txt
