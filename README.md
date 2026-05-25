# CIRI Analysis Pipeline — Local (No Docker)

## Output folder structure

Every step writes to its own **dated subfolder** under `Output/`, matching the catcheR convention:

```
Output/
  <YYMMDD>_step01_assignment_<sample>/
    csv/          tables
    plots/        PDF figures
    stats/        text summaries / run logs
    R_objects/    .RData objects (permanent archive)
    to_scratch/   files consumed by the NEXT step (auto-resolved)
  <YYMMDD>_step02_filter_<sample>/
    ...
```

`to_scratch/` acts as the handshake between steps. Each step automatically
finds the most recent `to_scratch/` from the previous step — you don't need
to pass file paths manually.

---

## Pipeline at a Glance

```
Step 00  Download Ensembl reference      (once)   → to_scratch/ensembl_protein_coding_genes.csv
Step 01  Perturbation assignment                   → to_scratch/annotation_data.csv
Step 02  Annotation & filtering                    → to_scratch/annotated_matrix.csv
Step 03  Load & preprocess (Monocle3, UMAP)        → to_scratch/processed_cds.RData
Step 04  Target validation (violins)               → plots/, csv/
Step 05  Cluster enrichment                        → plots/, csv/
Step 06  Subclustering & trajectory                → to_scratch/processed_cds_<group>.RData
                                                     to_scratch/pseudotime_<group>.csv
Step 07  Pseudotime statistics (KS test)           → plots/, csv/
Step 08  Signature scoring                         → plots/, csv/
```

Steps 04–05 and 07–08 do not write to `to_scratch/` because they are
terminal analyses — they produce final figures and tables only.

---

## Quickstart

### 1. Clone the repo and open in RStudio

```bash
git clone https://github.com/yourname/CIRI_pipeline.git
```

Open `CIRI_pipeline.Rproj` (or set Working Directory to the repo root in RStudio).

### 2. Copy `run_analysis.R` and fill in your paths

```bash
cp run_analysis_template.R run_analysis_AB011.R
```

Edit the `DATA_DIR`, `SAMPLE`, and step-specific parameters at the top, then
run each `system2()` block one at a time from the RStudio console.

### 3. Run from the RStudio Terminal

```bash
Rscript steps/step_01_perturbation_assignment.R \
  --data_dir /path/to/data \
  --matrix   filtered_feature_bc_matrix.h5 \
  --sample   AB011
```

Or from the R console:
```r
system2("Rscript", c("steps/step_01_perturbation_assignment.R",
  "--data_dir", "/path/to/data",
  "--matrix",   "filtered_feature_bc_matrix.h5",
  "--sample",   "AB011"))
```

---

## Step Reference

### Step 00 — Download Ensembl Gene Reference
```bash
Rscript steps/step_00_download_gene_ref.R --sample AB011
# output: Output/<date>_step00_geneRef_AB011/csv/ensembl_protein_coding_genes.csv
```

### Step 01 — Perturbation Assignment
```bash
Rscript steps/step_01_perturbation_assignment.R \
  --data_dir /path/to/data \
  --matrix   filtered_feature_bc_matrix.h5 \
  --sample   AB011 \
  --strategy 1          # 1=single guide, 2=dual guides
  # --threshold_a -1    # -1 = auto-detect via KDE valley
  # --threshold_i -1
```
Inspect: `Output/<date>_step01_assignment_AB011/plots/`

### Step 02 — Annotation & Filtering
```bash
Rscript steps/step_02_anno_filter.R \
  --data_dir /path/to/data \
  --matrix   filtered_feature_bc_matrix.h5 \
  --sample   AB011
  # --min_genes 250
  # --remove_mt TRUE   # explicitly removes ^MT- genes after protein-coding filter
  # --remove_rb TRUE
```
Inspect: `RiboMito_pre_filter.pdf`, `RiboMito_post_cell_filter.pdf`, `filter_summary.txt`

### Step 03 — Load & Preprocess
```bash
Rscript steps/step_03_load.R \
  --sample     AB011 \
  --resolution 5e-5    # increase → more clusters; decrease → fewer
```
Inspect: `umap_clusters.pdf`, `cluster_summary.txt`

### Step 04 — Target Validation
```bash
Rscript steps/step_04_target_validation.R \
  --sample    AB011 \
  --control_a "NTCa-NA"
```
Inspect: `validation_violin_<gene>.pdf`, `validation_summary.csv`

### Step 05 — Cluster Enrichment
```bash
Rscript steps/step_05_cluster_enrichment.R \
  --sample   AB011 \
  --clusters 5         # from umap_clusters.pdf; use "3-4" for multiple
  --control  "NTCa-NA"
```
Inspect: `heatmap_enrichment_Group_5.pdf`, `scatter_enrichment_Group_5.pdf`

### Step 06 — Subclustering & Trajectory
```bash
Rscript steps/step_06_subclusters_trajectory.R \
  --sample     AB011 \
  --clusters   5 \
  --root_gene  SOX2 \
  --group      muscle \
  --resolution 1e-3
```
Inspect: `UMAP_muscle_subclusters.pdf` → tune `--resolution`, then `UMAP_muscle_pseudotime.pdf`

### Step 07 — Pseudotime Statistics
```bash
Rscript steps/step_07_pseudotime_stats.R \
  --sample   AB011 \
  --group    muscle \
  --control  "NTCa-NA" \
  --min_cells 8
```
Inspect: `volcano_muscle.pdf`, `ks_results_muscle.csv`

### Step 08 — Signature Scoring
```bash
Rscript steps/step_08_signatures.R \
  --sample AB011 \
  --group  muscle
```
To add gene sets, edit the `SIGNATURES` list at the top of `step_08_signatures.R`.

---

## Input Files

| File | Description |
|------|-------------|
| `filtered_feature_bc_matrix.h5` | 10x Genomics output |
| `guides.csv` | Guide library: `feature, type (a/i), fixed (f/v)` — no header |

---

## R Packages

```r
install.packages(c(
  "dplyr","tidyr","stringr","purrr","ggplot2","quantmod","pracma",
  "hdf5r","Matrix","zoo","scales","Seurat","data.table",
  "viridis","ggrepel","gtools","colorspace","ggsignif","ggtext","biomaRt"
))
BiocManager::install(c("biomaRt","SingleCellExperiment","SummarizedExperiment"))
remotes::install_github("cole-trapnell-lab/monocle3")
```
