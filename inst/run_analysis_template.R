# =============================================================================
# CIRI Analysis Template
# =============================================================================
# Copy and rename this file per experiment:
#   cp run_analysis_template.R run_analysis_AB014_AB016.R
#
# Run each block ONE AT A TIME from the RStudio console.
# Inspect output after each step before promoting to scratch/.
#
# SETUP — choose one:
# =============================================================================

## Option A: installed from GitHub
# devtools::install_github("ebalmas/CIRIseq", ref = "V3")
library(CIRI)

## Option B: local development (source all R/ files)
# invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))

# =============================================================================
# Configuration — fill these in for your experiment
# =============================================================================
DATA_DIR    <- "/path/to/your/data"           # folder with H5, aggregation.csv
MATRIX      <- "filtered_feature_bc_matrix.h5"
AGGR_CSV    <- "aggregation.csv"              # from CellRanger aggr
GUIDES      <- "guides_harmonised.csv"        # produced by harmonise_guide_names.R
PC_RDATA    <- ""                             # scratch/scratch_protein_coding_genes.RData or ""
SAMPLE      <- "AB011"                        # used in all output folder names
OUTPUT_ROOT <- "Output"
SCRATCH_DIR <- "scratch"

# =============================================================================
# PRE-STEP — Harmonise guide names (run ONCE before Step 01)
# =============================================================================
# CellRanger collapses guide replicate names (e.g. ATF7IP_1A + ATF7IP_1B
# become ATF7IP_1 in the H5). This maps your guides.csv to the names
# actually present in the H5 so Step 01 can match them correctly.
#
# Inspect scratch/name_mapping.csv afterwards to verify all guides matched.
ciri_harmonise_guides(
  guides_path      = file.path(DATA_DIR, "scratch/guides_2.csv"),
  protospacer_path = file.path(DATA_DIR, "scratch/protospacer_calls_per_cell.csv"),
  out_path         = file.path(DATA_DIR, "scratch/guides_harmonised.csv")
)
# Inspect: scratch/name_mapping.csv   — check all guides matched correctly
# Then set GUIDES <- "scratch/guides_harmonised.csv" in the config above.

# =============================================================================
# STEP 00 — Download Ensembl gene reference (once per project)
# =============================================================================
ciri_step00_download_ref(
  output_root = OUTPUT_ROOT,
  sample      = SAMPLE
)
# Inspect: Output/<date>_step00_geneRef_<sample>/csv/ensembl_protein_coding_genes.csv
ciri_promote_scratch("step00_geneRef", sample = SAMPLE,
                     output_root = OUTPUT_ROOT, scratch_dir = SCRATCH_DIR)

# =============================================================================
# STEP 01 — Perturbation Assignment
# =============================================================================
# Uses the harmonised guides CSV and the H5 matrix.
# KDE threshold plots are saved to plots/ so you can verify the cut visually.
ciri_step01_assignment(
  data_dir    = DATA_DIR,
  matrix      = MATRIX,
  guides      = GUIDES,
  sample      = SAMPLE,
  output_root = OUTPUT_ROOT,
  strategy    = 2,       # 1 = single arm (CRISPRa OR CRISPRi)
                          # 2 = CIRI dual arm (CRISPRa AND CRISPRi)
  threshold_a = -1,      # -1 = auto-detect via KDE valley
  threshold_i = -1       # set manually if KDE plot looks wrong
)
# Inspect before promoting:
#   Output/<date>_step01_assignment_<sample>/plots/threshold_kde_CRISPRa.pdf
#   Output/<date>_step01_assignment_<sample>/plots/threshold_kde_CRISPRi.pdf
#   Output/<date>_step01_assignment_<sample>/plots/total_umixguide.pdf
#   Output/<date>_step01_assignment_<sample>/stats/assignment_summary.txt
#
# If thresholds look wrong, re-run with manual values, e.g.:
#   threshold_a = 140, threshold_i = 50
ciri_promote_scratch("step01_assignment", sample = SAMPLE,
                     output_root = OUTPUT_ROOT, scratch_dir = SCRATCH_DIR)

# =============================================================================
# STEP 02a — Annotate & QC  [NO FILTERING]
# =============================================================================
# Builds the Seurat object, joins aggregation.csv metadata + guide assignments,
# computes mito/ribo %, saves QC plots with suggested threshold dotted lines.
# Nothing is filtered here — just look at the plots and decide your real cuts.
ciri_step02_annotate(
  data_dir             = DATA_DIR,
  matrix               = "scratch/filtered_feature_bc_matrix.h5",
  aggr_csv             = "scratch/aggregation.csv",
  scratch_dir          = SCRATCH_DIR,
  output_root          = OUTPUT_ROOT,
  sample               = SAMPLE,
  protein_coding_rdata = PC_RDATA,
  # Suggested dotted lines (adjust to move lines on plots, not actual filters):
  suggest_mito_hi  = 15, suggest_mito_lo  = 1, suggest_ribo_lo  = 3,
  suggest_nGene_lo = 300, suggest_nGene_hi = 7000, suggest_nUMI_lo = 100
)
# ↑ Inspect before continuing:
#   Output/QC/<date>/ribomito/<date>_all_QC_final.pdf
#   Output/QC/<date>/ribomito/<date>_percent_MT_RIBO_QC.pdf
#   Output/QC/<date>/ribomito/<date>_gene_UMI_QC.pdf
#   Output/QC/<date>/csv/<date>_pre_QCcounts.csv
#
# Copy seurat_annotated.rds to scratch/ when happy:
file.copy("Output/QC/<date>/to_scratch/seurat_annotated.rds", SCRATCH_DIR)

# =============================================================================
# STEP 02b — Filter + Monocle3   (run up to 3 times)
# =============================================================================
# RUN 1 — apply your real filter thresholds:
ciri_step02_filter(
  scratch_dir  = SCRATCH_DIR,
  output_root  = OUTPUT_ROOT,
  sample       = SAMPLE,
  mito_hi  = 10, mito_lo = 0, ribo_lo = 1,
  nGene_lo = 300, nGene_hi = 7000, nUMI_lo = 100,
  num_dim  = 30
)
# ↑ Inspect: Output/monocle/<date>/umap/<date>_variance_knee_plot_dim30.png
#   → Curve should flatten before dim 30; if not, re-run with num_dim = 50
#
# RUN 2 (only if knee plot wrong — increase num_dim):
# ciri_step02_filter(scratch_dir=SCRATCH_DIR, output_root=OUTPUT_ROOT,
#   sample=SAMPLE, num_dim=50,
#   mito_hi=10, mito_lo=0, ribo_lo=1, nGene_lo=300, nGene_hi=7000, nUMI_lo=100)
#
# ↑ Compare 7 UMAP PDFs in Output/monocle/<date>/umap/:
#   _cds_6_ res=1e-5   (fewest) … _cds_1_ res=1e-2 (most clusters)
#
# RUN 3 — finalise chosen resolution:
# ciri_step02_filter(scratch_dir=SCRATCH_DIR, output_root=OUTPUT_ROOT,
#   sample=SAMPLE, chosen_cds="cds_3", num_dim=30,
#   mito_hi=10, mito_lo=0, ribo_lo=1, nGene_lo=300, nGene_hi=7000, nUMI_lo=100)
#
# Copy cds_final.rds to scratch/ when happy:
file.copy("Output/monocle/<date>/to_scratch/cds_final.rds", SCRATCH_DIR)

# =============================================================================
# STEP 03 — Load & Preprocess  (if not using step 02b's Monocle3 output)
# =============================================================================
# Only needed if you want to build the CDS from an annotated_matrix.csv
# rather than from the Seurat/H5 workflow above.
ciri_step03_load(
  sample       = SAMPLE,
  output_root  = OUTPUT_ROOT,
  scratch_dir  = SCRATCH_DIR,
  resolution   = 5e-5
)
ciri_promote_scratch("step03_load", sample = SAMPLE,
                     output_root = OUTPUT_ROOT, scratch_dir = SCRATCH_DIR)

# =============================================================================
# STEP 04 — Target Validation
# =============================================================================
ciri_step04_validation(
  sample      = SAMPLE,
  output_root = OUTPUT_ROOT,
  scratch_dir = SCRATCH_DIR,
  control_a   = "NTCa-NA",
  control_i   = "NTCi-NA"
)
# Inspect: Output/<date>_step04_validation_<sample>/plots/validation_violin_<gene>.pdf
#          Output/<date>_step04_validation_<sample>/csv/validation_summary.csv

# =============================================================================
# STEP 05 — Cluster Enrichment
# =============================================================================
TARGET_CLUSTER <- c("5")  # from UMAP; use c("3","4") for multiple

ciri_step05_enrichment(
  clusters    = TARGET_CLUSTER,
  control     = "NTCa-NA",
  sample      = SAMPLE,
  output_root = OUTPUT_ROOT,
  scratch_dir = SCRATCH_DIR,
  min_cells   = 10
)
# Inspect: Output/<date>_step05_enrichment_<sample>/plots/heatmap_enrichment_Group_5.pdf

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
  resolution  = 1e-3
)
# Inspect: Output/<date>_step06_trajectory_<sample>/plots/UMAP_muscle_subclusters.pdf
#          Output/<date>_step06_trajectory_<sample>/plots/UMAP_muscle_pseudotime.pdf
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
# Inspect: Output/<date>_step07_pseudotime_<sample>/plots/volcano_muscle.pdf
#          Output/<date>_step07_pseudotime_<sample>/csv/ks_results_muscle.csv

# =============================================================================
# STEP 08 — Signature Scoring
# =============================================================================
ciri_step08_signatures(
  group       = GROUP,
  sample      = SAMPLE,
  output_root = OUTPUT_ROOT,
  scratch_dir = SCRATCH_DIR
  # signatures = list(my_sig = c("MYOD1","MYOG","MYF5")) # add custom gene sets
)
# Inspect: Output/<date>_step08_signatures_<sample>/plots/umap_<sig>.pdf
#          Output/<date>_step08_signatures_<sample>/csv/signature_summary.csv
#          Output/<date>_step08_signatures_<sample>/stats/signature_coverage.txt
