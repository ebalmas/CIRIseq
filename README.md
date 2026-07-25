# CIRI — Single-Cell CRISPRa/i Screen Analysis Pipeline

End-to-end R package for CIRI screens: guide name harmonisation, perturbation
assignment with KDE auto-thresholding, Seurat v5 QC, Monocle3 preprocessing
at multiple clustering resolutions, cluster enrichment, pseudotime statistics,
and gene signature scoring.

Full dockerized version on Rstudio server is available from image hedgelab/rstudio-hedgelab:iPS2seq_CIRI_new_cellranger9

To run the Docker by command line:
```bash
docker run -d -p 8080:8787 --privileged=true --name container_name hedgelab/rstudio-hedgelab:iPS2seq_CIRI_new_cellranger9
docker exec -it container_name /bin/bash
```

To run the Docker in Rstudio server:
Note that USER and PASSWORD are the credential for Rstudio
```bash
docker run -d -itv /path/to/shared/folder:/scratch \

  --privileged=true \

  -p 8080:8787 \

  -e USER=rstudio \

  -e PASSWORD=your_password_here \

  --name container_name \

  hedgelab/rstudio-hedgelab:iPS2seq_CIRI_new_cellranger9

docker exec -idt container_name rstudio-server start
```
go on browser http://localhost:8080/ and add the credentials


## Below are the instructions to use outside the prepared Docker:
---

## Installation

CIRI itself installs in seconds. Heavy bioinformatics packages (`Seurat`,
`monocle3`, `hdf5r`) are **optional** — each step checks for what it needs
at runtime and gives a clear install instruction if something is missing.

### Step 1 — Install CIRI (requires R ≥ 4.1)

```r
install.packages("devtools")
library(devtools)
devtools::install_github("ebalmas/CIRIseq", ref = "V4")
library(CIRI)
```

### Step 2 — Install optional dependencies when ready

```r
# See what is installed and what is missing
check_dependencies()

# Install everything at once
install_dependencies()

# Or install only what you need for specific steps
install_dependencies(steps = "00")      # biomaRt  (Ensembl reference download)
install_dependencies(steps = "01")      # hdf5r    (reading .h5 files)
install_dependencies(steps = "02")      # hdf5r + Seurat + SeuratObject + patchwork
install_dependencies(steps = "03-08")   # monocle3 + igraph
```

`check_dependencies()` output:

```
CIRI dependency status:
--------------------------------------------------
  [OK]    Step 00 — Ensembl reference         ready
  [MISS]  Step 01 — Guide assignment          missing: hdf5r
  [MISS]  Step 02 — Seurat QC                 missing: Seurat, SeuratObject, patchwork
  [MISS]  Steps 03-08 — Monocle3              missing: monocle3, igraph
--------------------------------------------------
```

### Package versions

| Package | Minimum version | Notes |
|---|---|---|
| R | ≥ 4.1.0 | Required for native pipe `\|>` |
| Seurat | ≥ 5.0.0 | v5 changed `GetAssayData(slot=)` → `layer=` |
| SeuratObject | ≥ 5.0.0 | Ships separately from Seurat since v5 |
| hdf5r | ≥ 1.3.8 | Requires HDF5 system library (see below) |
| ggplot2 | ≥ 3.4.0 | Required for `after_stat()` |
| dplyr | ≥ 1.1.0 | Required for `case_when()` improvements |
| monocle3 | GitHub only | Install via `install_dependencies()` |

### hdf5r system library (macOS / Linux only)

Install before `install_dependencies(steps = "01")`:

```bash
brew install hdf5                      # macOS (Homebrew)
sudo apt-get install libhdf5-dev       # Ubuntu / Debian
sudo dnf install hdf5-devel            # Fedora / RHEL
```

### Seurat v5 note

Seurat v5 (released 2023) broke the `slot =` argument in `GetAssayData()`.
The CIRI package uses `SeuratObject::GetAssayData(assay = "RNA", layer = "counts")`
throughout. If you see:

```
Error: The `slot` argument of `GetAssayData()` was deprecated in SeuratObject 5.0.0
```

update the CIRI package: `devtools::install_github("ebalmas/CIRIseq", ref = "V4", force = TRUE)`

---

## Pipeline overview

```
Pre-step  ciri_harmonise_guides()    fix guide name mismatches vs H5
Step 00   ciri_step00_download_ref() download Ensembl protein-coding genes (once)
Step 01   ciri_step01_assignment()   assign CRISPRa/i guides to cells
Step 02a  ciri_step02_annotate()     Seurat object + QC plots (NO filtering)
Step 02b  ciri_step02_filter()       filter cells + Monocle3 CDS (3 runs)
Step 03   ciri_step03_load()         alternative CSV-based Monocle3 entry point
Step 04   ciri_step04_validation()   target knockdown/activation validation
Step 05   ciri_step05_enrichment()   cluster enrichment analysis
Step 06   ciri_step06_trajectory()   subclustering + pseudotime trajectory
Step 07   ciri_step07_pseudotime()   KS test on pseudotime distributions
Step 08   ciri_step08_signatures()   gene set signature scoring
```

Steps 02a and 02b are available **both as package functions** and as
**standalone scripts** (`ciri_copy_scripts()`) for terminal use.

---

## Folder structure

```
my_analysis/
├── scratch/                                ← staging area YOU control
│   ├── filtered_feature_bc_matrix.h5       from CellRanger aggr
│   ├── aggregation.csv                     input to CellRanger aggr
│   ├── protospacer_calls_per_cell.csv       from CellRanger aggr
│   ├── guides_2.csv                         your original guide library CSV
│   ├── guides_harmonised.csv               produced by ciri_harmonise_guides()
│   └── scratch_protein_coding_genes.RData  optional — protein-coding gene filter
├── Output/                                 ← dated run folders (auto-created)
│   ├── QC/<YYMMDD>/
│   │   ├── csv/   ribomito/   filtering/   R_objects/
│   │   └── to_scratch/
│   ├── monocle/<YYMMDD>/
│   │   ├── umap/   R_objects/
│   │   └── to_scratch/
│   └── <YYMMDD>_step01_assignment_<sample>/
│       ├── csv/   plots/   stats/   R_objects/
│       └── to_scratch/
├── harmonise_guide_names.R                 copied by ciri_copy_scripts()
├── ciri_step02_annotate.R                  copied by ciri_copy_scripts()
├── ciri_step02_filter.R                    copied by ciri_copy_scripts()
└── run_analysis_AB011.R                    your copy of the template
```

Add `Output/`, `scratch/*.h5`, and `scratch/*.rds` to `.gitignore`.
Commit your `run_analysis_<experiment>.R` — it is your lab notebook.

---

## Quick start

```r
library(CIRI)

# 1. Copy scripts and analysis template to your working directory
ciri_copy_scripts()
file.copy(system.file("run_analysis_template.R", package = "CIRI"),
          "run_analysis_AB011.R")

# 2. Open run_analysis_AB011.R and fill in DATA_DIR, SAMPLE, etc.
# 3. Run each block one at a time from the RStudio console.
```

---

## The scratch/ checkpoint system

Nothing moves to `scratch/` automatically. After each step you inspect the
output and decide when to promote:

```r
ciri_step01_assignment(data_dir = "/data", matrix = "matrix.h5", sample = "AB011")

# Inspect Output/<date>_step01_assignment_AB011/plots/ and stats/
# Happy with the result? Promote to scratch/:
ciri_promote_scratch("step01_assignment", sample = "AB011")

# Next step reads from scratch/ automatically
```

| Function | Description |
|---|---|
| `ciri_promote_scratch()` | Copy `to_scratch/` → `scratch/` |
| `ciri_list_scratch()` | List what is ready in a step's `to_scratch/` |
| `ciri_scratch_status()` | Show current `scratch/` contents |
| `ciri_clear_scratch()` | Clear `scratch/` (requires `confirm = TRUE`) |
| `ciri_copy_scripts()` | Copy standalone scripts to working directory |

---

## Pre-step — Guide name harmonisation

CellRanger collapses guide replicate names in the H5 (e.g. `ATF7IP_1A` +
`ATF7IP_1B` → `ATF7IP_1`). Without harmonisation, Step 01 cannot match the
guides CSV to the H5 features and returns 0 assigned cells.

```r
ciri_harmonise_guides(
  guides_path      = "scratch/guides_2.csv",
  protospacer_path = "scratch/protospacer_calls_per_cell.csv",
  out_path         = "scratch/guides_harmonised.csv"
)
```

Inspect `scratch/name_mapping.csv` to verify all guides matched. Then use
`guides_harmonised.csv` in Step 01.

---

## Step 01 — Perturbation assignment

```r
ciri_step01_assignment(
  data_dir    = "/path/to/data",
  matrix      = "scratch/filtered_feature_bc_matrix.h5",
  guides      = "scratch/guides_harmonised.csv",
  sample      = "AB011",
  strategy    = 2,       # 1 = single arm, 2 = CIRI dual arm
  threshold_a = -1,      # -1 = auto-detect via KDE valley
  threshold_i = -1
)
```

**KDE auto-thresholding:** when `threshold = -1`, the pipeline finds the
first valley in the kernel density estimate of per-cell fixed-guide UMIs —
the dip between the noise peak (cells without the guide) and the signal peak
(cells with it). Two diagnostic PDFs are saved to `plots/` so you can verify:

- `plots/threshold_kde_CRISPRa.pdf` — red dashed line = chosen threshold
- `plots/threshold_kde_CRISPRi.pdf`

If the line is in the wrong place, override with `threshold_a = 140`.

---

## Step 02a — Annotate & QC

```r
ciri_step02_annotate(
  data_dir             = "/path/to/data",
  matrix               = "scratch/filtered_feature_bc_matrix.h5",
  aggr_csv             = "scratch/aggregation.csv",
  scratch_dir          = "scratch",
  sample               = "AB011",
  protein_coding_rdata = "",  # optional protein-coding gene RData
  # Suggested dotted lines on plots (not actual filters):
  suggest_mito_hi  = 15, suggest_mito_lo  = 1, suggest_ribo_lo  = 3,
  suggest_nGene_lo = 300, suggest_nGene_hi = 7000, suggest_nUMI_lo = 100
)
```

Builds the Seurat object, joins aggregation metadata and guide assignments,
computes `percent.mt` and `percent.ribo`, saves 8 QC PDFs. **No cells are
removed.** Threshold values are dotted lines only — inspect the plots and
decide your real cuts before running Step 02b.

Or use the standalone script from the terminal:

```bash
Rscript ciri_step02_annotate.R \
  --data_dir /path/to/data --sample AB011 \
  --mito_hi 15 --nGene_lo 300 --nUMI_lo 100
```

---

## Step 02b — Filter + Monocle3

Run up to three times:

**Run 1** — apply thresholds, build CDS at all 7 resolutions:

```r
ciri_step02_filter(
  scratch_dir = "scratch", sample = "AB011",
  mito_hi = 10, mito_lo = 0, ribo_lo = 1,
  nGene_lo = 300, nGene_hi = 7000, nUMI_lo = 100,
  num_dim = 30
)
```

→ Inspect `Output/monocle/<date>/umap/<date>_variance_knee_plot_dim30.png`
If the variance curve has not flattened by dim 30, re-run with `num_dim = 50`.

**Run 2** (if knee plot wrong):

```r
ciri_step02_filter(...same args..., num_dim = 50)
```

→ Compare 7 UMAP PDFs in `Output/monocle/<date>/umap/`:

| File suffix | Resolution | Expected |
|---|---|---|
| `_cds_6_` | `1e-5` | fewest clusters |
| `_cds_7_` | `2.5e-5` | very coarse |
| `_cds_2_` | `1e-4` | coarse |
| `_cds_4_` | `2e-4` | moderate |
| `_cds_3_` | `2e-4, k=15` | moderate (default) |
| `_cds_5_` | `5e-4` | fine |
| `_cds_1_` | `1e-2` | finest — usually too many |

**Run 3** — finalise chosen resolution:

```r
ciri_step02_filter(...same args..., chosen_cds = "cds_3", num_dim = 30)
```

→ Saves `cds_final.rds` to `Output/monocle/<date>/to_scratch/`.

Or use the standalone script equivalently:

```bash
# Run 1
Rscript ciri_step02_filter.R --sample AB011 --mito_hi 10 --nGene_lo 300
# Run 3
Rscript ciri_step02_filter.R --sample AB011 --chosen_cds cds_3 --num_dim 30 \
  --mito_hi 10 --mito_lo 0 --ribo_lo 1 --nGene_lo 300 --nGene_hi 7000 --nUMI_lo 100
```

---

## Parameter reference

### `ciri_harmonise_guides()`

| Parameter | Description |
|---|---|
| `guides_path` | Original guides CSV (no header: feature, type, fixed) |
| `protospacer_path` | `protospacer_calls_per_cell.csv` from CellRanger |
| `out_path` | Output path. Default: `guides_harmonised.csv` next to guides |

### `ciri_step01_assignment()`

| Parameter | Default | Description |
|---|---|---|
| `data_dir` | required | Folder with H5 and guides CSV |
| `matrix` | required | H5 filename |
| `guides` | `"guides.csv"` | Use `guides_harmonised.csv` |
| `sample` | `"CIRI"` | Experiment name |
| `strategy` | `1` | `1` = single arm; `2` = CIRI dual arm |
| `threshold_a` | `-1` | CRISPRa UMI threshold; `-1` = auto KDE |
| `threshold_i` | `-1` | CRISPRi UMI threshold; `-1` = auto KDE |

### `ciri_step02_annotate()`

| Parameter | Default | Description |
|---|---|---|
| `data_dir` | required | Folder with H5 and aggregation.csv |
| `matrix` | `"scratch/filtered_feature_bc_matrix.h5"` | H5 filename |
| `aggr_csv` | `"scratch/aggregation.csv"` | CellRanger aggr input CSV |
| `scratch_dir` | `"scratch"` | Folder with annotation_data.csv |
| `protein_coding_rdata` | `""` | Optional gene filter RData |
| `suggest_mito_hi/lo` | `15` / `1` | Dotted lines on plots only |
| `suggest_ribo_lo` | `3` | Dotted line on plots only |
| `suggest_nGene_lo/hi` | `300` / `7000` | Dotted lines on plots only |
| `suggest_nUMI_lo` | `100` | Dotted line on plots only |

### `ciri_step02_filter()`

| Parameter | Default | Description |
|---|---|---|
| `scratch_dir` | `"scratch"` | Folder with `seurat_annotated.rds` |
| `mito_lo` / `mito_hi` | `0` / `10` | `percent.mt` bounds |
| `ribo_lo` | `1` | `percent.ribo` lower bound |
| `nGene_lo` / `nGene_hi` | `300` / `7000` | `nFeature_RNA` bounds |
| `nUMI_lo` | `100` | `nCount_RNA` lower bound |
| `num_dim` | `30` | PCA dimensions for `preprocess_cds()` |
| `chosen_cds` | `"cds_3"` | Which resolution to save as final CDS |

### `ciri_step06_trajectory()`

| Parameter | Default | Description |
|---|---|---|
| `clusters` | required | Cluster ID(s), e.g. `"5"` or `c("3","4")` |
| `root_gene` | required | Gene highest at trajectory start |
| `group` | required | Lineage name, e.g. `"muscle"` |
| `resolution` | `1e-3` | Sub-clustering resolution |
| `n_dims` | `50` | PCA dims for subset |

---

## Guides CSV format

No header. Three columns: `feature`, `type` (`a`/`i`), `fixed` (`f`/`v`).

```
MYOD_1,a,f      ← fixed CRISPRa barcoding guide
MYOD_2,a,f
NANOG,i,f       ← fixed CRISPRi barcoding guide
OCT4,i,f
SOX2,i,f
BAF60C_1,a,v    ← variable CRISPRa perturbation guide
BAF60C_2,a,v
CTCF_1,i,v      ← variable CRISPRi perturbation guide
NTCa,a,v        ← non-targeting control CRISPRa
NTCi,i,v        ← non-targeting control CRISPRi
```

**Always use `guides_harmonised.csv`** (from `ciri_harmonise_guides()`) with
Step 01. The names in this file match exactly what is in the H5.

---

## Key fixes vs original pipeline

| Issue | Fix |
|---|---|
| Guide name mismatch (H5 uses collapsed names, guides CSV has A/B suffixes) | `ciri_harmonise_guides()` pre-step |
| Step 01 returned 0 assigned cells | Fixed by guide name harmonisation |
| `GetAssayData(slot=)` defunct in Seurat ≥ 5.0.0 | Updated to `SeuratObject::GetAssayData(layer=)` |
| `Remotes:` field forced monocle3 install on `devtools::install_github()` | Removed from DESCRIPTION |
| NAMESPACE multi-line `importFrom` caused silent install failure | One `importFrom` per line |
| MT-encoded genes survive protein-coding biomaRt filter | `^MT-` genes explicitly removed after filter |
| `subset()` deprecated in Seurat v5 | Uses `base::subset()` on Seurat object |
