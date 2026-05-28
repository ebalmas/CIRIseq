# CIRI Pipeline — Step-by-Step Protocol

## What you need before starting

### Software
- R ≥ 4.1.0
- RStudio (recommended)
- HDF5 system library (for reading `.h5` files):
  ```bash
  brew install hdf5              # macOS
  sudo apt-get install libhdf5-dev   # Ubuntu/Debian
  sudo dnf install hdf5-devel        # Fedora/RHEL
  ```

### Files from CellRanger aggr output
All files go in your `scratch/` folder:

| File | Where to get it |
|---|---|
| `filtered_feature_bc_matrix.h5` | `<aggr_run>/outs/` |
| `aggregation.csv` | the CSV you passed to `cellranger aggr` |
| `protospacer_calls_per_cell.csv` | `<aggr_run>/outs/` |
| `guides_2.csv` | your guide library annotation (from lab) |
| `scratch_protein_coding_genes.RData` | optional — from previous runs |

---

## Folder setup

Create this structure in your project directory:

```
my_analysis/
└── scratch/
    ├── filtered_feature_bc_matrix.h5
    ├── aggregation.csv
    ├── protospacer_calls_per_cell.csv
    └── guides_2.csv
```

Open RStudio and set the working directory to `my_analysis/`:
```r
setwd("/path/to/my_analysis")
```

---

## Install the package

```r
install.packages("devtools")
devtools::install_github("ebalmas/CIRIseq", ref = "V3")
library(CIRI)

# Install bioinformatics dependencies
install_dependencies()   # installs everything; or use steps = "01" etc.

# Verify
check_dependencies()
```

---

## Copy scripts and create your analysis file

```r
library(CIRI)

# Copy the 3 standalone scripts and the analysis template
ciri_copy_scripts()   # creates harmonise_guide_names.R, ciri_step02_annotate.R,
                      #         ciri_step02_filter.R  in your working directory
file.copy(
  system.file("run_analysis_template.R", package = "CIRI"),
  "run_analysis_AB011.R"
)
```

Open `run_analysis_AB011.R` and fill in the configuration block at the top:
```r
DATA_DIR    <- "/path/to/my_analysis"
MATRIX      <- "filtered_feature_bc_matrix.h5"
AGGR_CSV    <- "aggregation.csv"
GUIDES      <- "guides_harmonised.csv"   # will exist after the pre-step
SAMPLE      <- "AB011"
OUTPUT_ROOT <- "Output"
SCRATCH_DIR <- "scratch"
```

---

## PRE-STEP — Harmonise guide names

CellRanger collapses guide replicate names in the H5 (e.g. `ATF7IP_1A` +
`ATF7IP_1B` → `ATF7IP_1`). This step produces a corrected guides CSV that
matches the H5 exactly.

**Run from the RStudio Terminal** (not the console):

```bash
Rscript harmonise_guide_names.R \
  --guides      scratch/guides_2.csv \
  --protospacer scratch/protospacer_calls_per_cell.csv \
  --out         scratch/guides_harmonised.csv
```

**Check the output:**
- `scratch/guides_harmonised.csv` — the corrected guide CSV
- `scratch/name_mapping.csv` — shows how each original name was matched

Look for any `UNMATCHED` warnings — these guides were not detected in the data
and will be excluded. This is expected for guides absent from the experiment.

---

## STEP 00 — Download Ensembl gene reference

Run **once per project** (needs internet). Saves a local CSV of human
protein-coding genes used to filter the expression matrix.

```r
ciri_step00_download_ref(output_root = OUTPUT_ROOT, sample = SAMPLE)
```

**Check:** `Output/<date>_step00_geneRef_<sample>/csv/ensembl_protein_coding_genes.csv`

```r
ciri_promote_scratch("step00_geneRef", sample = SAMPLE,
                     output_root = OUTPUT_ROOT, scratch_dir = SCRATCH_DIR)
```

---

## STEP 01 — Perturbation Assignment

Assigns CRISPRa and CRISPRi guide identities to each cell.

```r
ciri_step01_assignment(
  data_dir    = DATA_DIR,
  matrix      = "scratch/filtered_feature_bc_matrix.h5",
  guides      = "scratch/guides_harmonised.csv",
  sample      = SAMPLE,
  output_root = OUTPUT_ROOT,
  strategy    = 2,        # 2 = CIRI dual arm (CRISPRa AND CRISPRi required)
  threshold_a = -1,       # auto-detect
  threshold_i = -1
)
```

**Check these plots** before promoting:

```
Output/<date>_step01_assignment_<sample>/
  plots/
    threshold_kde_CRISPRa.pdf   ← red dashed line = chosen threshold
    threshold_kde_CRISPRi.pdf     should sit in the valley between noise and signal
    total_umixguide.pdf
    fixed_vs_variable_scatter.pdf
  stats/
    assignment_summary.txt       ← how many cells got each arm assigned
```

If the KDE threshold looks wrong, re-run with manual values:
```r
ciri_step01_assignment(..., threshold_a = 140, threshold_i = 50)
```

When happy:
```r
ciri_promote_scratch("step01_assignment", sample = SAMPLE,
                     output_root = OUTPUT_ROOT, scratch_dir = SCRATCH_DIR)
```

---

## STEP 02a — Annotate & QC plots

Builds the Seurat object and produces QC plots. **No cells are removed.**
The threshold values are shown as dotted lines so you can decide your real cuts.

**Run from the RStudio Terminal:**

```bash
Rscript ciri_step02_annotate.R \
  --data_dir /path/to/my_analysis \
  --sample   AB011 \
  --mito_hi  15  --mito_lo 1  --ribo_lo 3 \
  --nGene_lo 300 --nGene_hi 7000 --nUMI_lo 100
```

Or from RStudio: edit `INTERACTIVE_PARAMS` at the top of `ciri_step02_annotate.R`
and `source("ciri_step02_annotate.R")`.

**Check these plots:**

```
Output/QC/<date>/ribomito/
  <date>_all_QC_final.pdf          ← 3-panel: mito, genes, ribo scatter
  <date>_percent_MT_RIBO_QC.pdf    ← ribo vs mito scatter by sample
  <date>_gene_UMI_QC.pdf           ← nUMI vs nGene coloured by mito %
  <date>_feature_QC.pdf            ← violin plots per sample
Output/QC/<date>/csv/
  <date>_pre_QCcounts.csv          ← cells per sample before filtering
```

Decide your filter thresholds from these plots, then:
```r
file.copy("Output/QC/<date>/to_scratch/seurat_annotated.rds", "scratch/")
```

---

## STEP 02b — Filter + Monocle3

### Run 1 — apply thresholds, generate all UMAPs

```bash
Rscript ciri_step02_filter.R \
  --sample   AB011 \
  --mito_hi  10  --mito_lo 0  --ribo_lo 1 \
  --nGene_lo 300 --nGene_hi 7000 --nUMI_lo 100
```

**Check the knee plot:**
```
Output/monocle/<date>/umap/<date>_variance_knee_plot_dim30.png
```
The curve should flatten (level off) before dim 30. If it's still dropping
steeply at dim 30, re-run with more dimensions.

### Run 2 (only if knee plot says more dims needed)

```bash
Rscript ciri_step02_filter.R \
  --sample AB011 --num_dim 50 \
  --mito_hi 10 --mito_lo 0 --ribo_lo 1 \
  --nGene_lo 300 --nGene_hi 7000 --nUMI_lo 100
```

**Check the 7 UMAP resolutions** in `Output/monocle/<date>/umap/`:

| File suffix | Resolution | Expected clusters |
|---|---|---|
| `_cds_6_` | `1e-5` | fewest — very coarse |
| `_cds_7_` | `2.5e-5` | coarse |
| `_cds_2_` | `1e-4` | moderate |
| `_cds_4_` | `2e-4` | moderate |
| `_cds_3_` | `2e-4 k=15` | moderate (default) |
| `_cds_5_` | `5e-4` | fine |
| `_cds_1_` | `1e-2` | finest — usually too many |

Pick the resolution that best separates your biology.

### Run 3 — finalise

```bash
Rscript ciri_step02_filter.R \
  --sample      AB011 \
  --chosen_cds  cds_3 \
  --num_dim     30 \
  --mito_hi     10 --mito_lo 0 --ribo_lo 1 \
  --nGene_lo    300 --nGene_hi 7000 --nUMI_lo 100
```

Then:
```r
file.copy("Output/monocle/<date>/to_scratch/cds_final.rds", "scratch/")
```

---

## STEP 04 — Target Validation

Checks knockdown/activation efficiency for each targeted gene.

```r
ciri_step04_validation(
  sample      = SAMPLE,
  output_root = OUTPUT_ROOT,
  scratch_dir = SCRATCH_DIR,
  control_a   = "NTCa",    # match exactly to your NTC guide name
  control_i   = "NTCi"
)
```

**Check:** `Output/<date>_step04_validation_<sample>/plots/validation_violin_<gene>.pdf`

---

## STEP 05 — Cluster Enrichment

What fraction of each perturbation ends up in a target cluster?

```r
TARGET_CLUSTER <- c("5")   # from the UMAP — change to your target cluster(s)

ciri_step05_enrichment(
  clusters    = TARGET_CLUSTER,
  control     = "NTCa",
  sample      = SAMPLE,
  output_root = OUTPUT_ROOT,
  scratch_dir = SCRATCH_DIR,
  min_cells   = 10
)
```

**Check:** `Output/<date>_step05_enrichment_<sample>/plots/heatmap_enrichment_Group_5.pdf`

---

## STEP 06 — Subclustering & Trajectory

Subsets the target cluster, re-clusters, and learns a pseudotime trajectory.

```r
ciri_step06_trajectory(
  clusters    = TARGET_CLUSTER,
  root_gene   = "SOX2",     # gene highest at the start of differentiation
  group       = "muscle",
  sample      = SAMPLE,
  output_root = OUTPUT_ROOT,
  scratch_dir = SCRATCH_DIR,
  resolution  = 1e-3
)
```

**Check before promoting:**
```
Output/<date>_step06_trajectory_<sample>/plots/
  UMAP_muscle_subclusters.pdf   ← tune --resolution if needed
  UMAP_muscle_pseudotime.pdf    ← check root placement
```

```r
ciri_promote_scratch("step06_trajectory", sample = SAMPLE,
                     output_root = OUTPUT_ROOT, scratch_dir = SCRATCH_DIR)
```

---

## STEP 07 — Pseudotime Statistics

KS test: which perturbations shift pseudotime vs the control?

```r
ciri_step07_pseudotime(
  group       = "muscle",
  control     = "NTCa-NA",   # gene_comb format: <gene_a>-<gene_i>
  sample      = SAMPLE,
  output_root = OUTPUT_ROOT,
  scratch_dir = SCRATCH_DIR,
  min_cells   = 8,
  ecdf_top_n  = 10
)
```

**Check:**
```
Output/<date>_step07_pseudotime_<sample>/
  plots/volcano_muscle.pdf         ← KS stat vs -log10(p.adj)
  csv/ks_results_muscle.csv        ← full results table
```

---

## STEP 08 — Signature Scoring

Scores cells against gene set signatures (cell cycle, myogenic, pluripotency, etc.).

```r
ciri_step08_signatures(
  group       = "muscle",
  sample      = SAMPLE,
  output_root = OUTPUT_ROOT,
  scratch_dir = SCRATCH_DIR
)
```

To add your own gene sets:
```r
ciri_step08_signatures(
  ...,
  signatures = list(
    my_sig_1 = c("GENE1", "GENE2", "GENE3"),
    my_sig_2 = c("GENEA", "GENEB")
  )
)
```

**Check:**
```
Output/<date>_step08_signatures_<sample>/
  stats/signature_coverage.txt     ← how many genes from each set were found
  plots/umap_<sig>.pdf             ← one UMAP per signature
  csv/signature_summary.csv        ← mean score per perturbation per signature
```

---

## Troubleshooting

| Problem | Solution |
|---|---|
| `Package 'monocle3' not installed` | `install_dependencies(steps = "03-08")` |
| `hdf5r` install fails | Install HDF5 system library first (see top of protocol) |
| Step 01 assigns 0 cells | Run `harmonise_guide_names.R` first; use `guides_harmonised.csv` |
| KDE threshold wrong | Pass `threshold_a = X` manually; inspect `threshold_kde_CRISPRa.pdf` |
| `GetAssayData(slot=)` error | Update CIRI: `devtools::install_github("ebalmas/CIRIseq", ref = "V3", force = TRUE)` |
| Knee plot curve still dropping at dim 30 | Re-run `ciri_step02_filter.R` with `--num_dim 50` |
| Too many / too few clusters | Change `--chosen_cds` in step 02b Run 3 |
