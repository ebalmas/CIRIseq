# CIRI — Single-Cell CRISPRa/i Screen Analysis Pipeline

End-to-end R package for CIRI screens: guide name harmonisation, perturbation
assignment, Seurat QC, Monocle3 preprocessing, cluster enrichment, pseudotime
statistics, and gene signature scoring.

---

## Installation

CIRI installs in seconds — heavy packages (`Seurat`, `monocle3`, `hdf5r`) are
optional and only loaded when the step that needs them is called.

```r
install.packages("devtools")
devtools::install_github("ebalmas/CIRIseq", ref = "V3")
library(CIRI)

# Check what optional dependencies are installed
check_dependencies()

# Install everything at once
install_dependencies()

# Or install per step
install_dependencies(steps = "01")      # hdf5r  (H5 reading)
install_dependencies(steps = "02")      # hdf5r + Seurat + data.table
install_dependencies(steps = "03-08")   # monocle3 + igraph
install_dependencies(steps = "00")      # biomaRt (Ensembl reference)
```

### hdf5r system library (required before installing hdf5r)

```bash
brew install hdf5           # macOS
sudo apt-get install libhdf5-dev   # Ubuntu/Debian
sudo dnf install hdf5-devel        # Fedora/RHEL
```

---

## Pipeline overview

```
Pre-step  harmonise_guide_names.R    fix guide name mismatches vs H5
Step 00   ciri_step00_download_ref() download Ensembl protein-coding genes
Step 01   ciri_step01_assignment()   assign CRISPRa/i guides to cells
Step 02a  ciri_step02_annotate.R     Seurat object + QC plots (standalone)
Step 02b  ciri_step02_filter.R       filter cells + Monocle3 CDS (standalone)
Step 03   ciri_step03_load()         (alternative to 02b for CSV-based input)
Step 04   ciri_step04_validation()   target knockdown/activation validation
Step 05   ciri_step05_enrichment()   cluster enrichment analysis
Step 06   ciri_step06_trajectory()   subclustering + pseudotime trajectory
Step 07   ciri_step07_pseudotime()   KS test on pseudotime distributions
Step 08   ciri_step08_signatures()   gene set signature scoring
```

Steps 02a and 02b are **standalone scripts** (not package functions) because
the Seurat QC step is intentionally kept separate — you can substitute your
own QC pipeline before the Monocle3 steps.

---

## Folder structure

```
my_analysis/
├── scratch/                          ← staging area YOU control
│   ├── filtered_feature_bc_matrix.h5
│   ├── aggregation.csv
│   ├── guides_2.csv
│   ├── guides_harmonised.csv         ← produced by pre-step
│   ├── protospacer_calls_per_cell.csv
│   └── scratch_protein_coding_genes.RData   (optional)
├── Output/                           ← dated run folders (auto-created)
│   ├── QC/
│   │   └── <YYMMDD>/
│   │       ├── csv/   ribomito/   filtering/   R_objects/
│   │       └── to_scratch/
│   ├── monocle/
│   │   └── <YYMMDD>/
│   │       ├── umap/   R_objects/   dotplot/
│   │       └── to_scratch/
│   └── <YYMMDD>_step01_assignment_<sample>/
│       ├── csv/   plots/   stats/   R_objects/
│       └── to_scratch/
├── harmonise_guide_names.R           ← copied by ciri_copy_scripts()
├── ciri_step02_annotate.R            ← copied by ciri_copy_scripts()
├── ciri_step02_filter.R              ← copied by ciri_copy_scripts()
└── run_analysis_AB011.R              ← your copy of the template
```

Add `Output/` and `scratch/*.h5` to `.gitignore`. Commit your
`run_analysis_<experiment>.R` — it is your lab notebook.

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
# Happy with the result? Promote it:
ciri_promote_scratch("step01_assignment", sample = "AB011")

# Next step reads from scratch/ automatically
```

Helper functions:

| Function | Description |
|---|---|
| `ciri_promote_scratch()` | Copy `to_scratch/` → `scratch/` |
| `ciri_list_scratch()` | List what's ready in a step's `to_scratch/` |
| `ciri_scratch_status()` | Show current `scratch/` contents |
| `ciri_clear_scratch()` | Clear `scratch/` (requires `confirm = TRUE`) |
| `ciri_copy_scripts()` | Copy standalone QC scripts to working directory |

---

## Step 01 — guide name harmonisation

CellRanger collapses guide replicate names: `ATF7IP_1A` + `ATF7IP_1B`
become `ATF7IP_1` in the H5. Run `harmonise_guide_names.R` once to
produce a `guides_harmonised.csv` that matches the H5 exactly.

```bash
Rscript harmonise_guide_names.R \
  --guides      scratch/guides_2.csv \
  --protospacer scratch/protospacer_calls_per_cell.csv \
  --out         scratch/guides_harmonised.csv
```

Inspect `scratch/name_mapping.csv` to verify all guides matched.

### Auto-thresholding

When `threshold_a = -1` (default), the threshold is auto-detected using the
first valley in the KDE of per-cell fixed-guide UMI sums — the dip between
the noise peak (cells that didn't receive the fixed guide) and the signal
peak (cells that did).

Two diagnostic plots are saved to `plots/`:
- `threshold_kde_CRISPRa.pdf`
- `threshold_kde_CRISPRi.pdf`

If the red dashed line lands in the wrong place, pass the threshold manually:

```r
ciri_step01_assignment(..., threshold_a = 140, threshold_i = 50)
```

---

## Step 02 — QC and filtering (standalone scripts)

### 02a — Annotate & QC plots

```bash
Rscript ciri_step02_annotate.R \
  --data_dir /path/to/data \
  --sample   AB011 \
  --mito_hi  15  --mito_lo 1  --ribo_lo 3 \
  --nGene_lo 300 --nGene_hi 7000 --nUMI_lo 100
```

Threshold values are shown as **dotted lines only** — no cells are removed.
Inspect the PDFs in `Output/QC/<date>/ribomito/` then decide your real cuts.

### 02b — Filter + Monocle3 (3 runs)

**Run 1** — apply thresholds, build CDS at 7 resolutions:
```bash
Rscript ciri_step02_filter.R \
  --sample   AB011 \
  --mito_hi  10  --mito_lo 0  --ribo_lo 1 \
  --nGene_lo 300 --nGene_hi 7000 --nUMI_lo 100
```
→ Inspect knee plot (`_variance_knee_plot_dim30.png`)

**Run 2** — if knee plot says more dims needed:
```bash
Rscript ciri_step02_filter.R --sample AB011 --num_dim 50 [+ same thresholds]
```
→ Inspect 7 UMAP PDFs (`_cds_6_` = fewest clusters … `_cds_1_` = most)

**Run 3** — finalise chosen resolution:
```bash
Rscript ciri_step02_filter.R \
  --sample AB011 --chosen_cds cds_3 --num_dim 30 \
  --mito_hi 10 --mito_lo 0 --ribo_lo 1 \
  --nGene_lo 300 --nGene_hi 7000 --nUMI_lo 100
```

---

## Parameter reference

### `ciri_step01_assignment()`

| Parameter | Default | Description |
|---|---|---|
| `data_dir` | required | Folder with H5 matrix and guides CSV |
| `matrix` | required | H5 filename |
| `guides` | `"guides.csv"` | Guide annotation CSV — use `guides_harmonised.csv` |
| `sample` | `"CIRI"` | Experiment name |
| `strategy` | `1` | `1` = single arm; `2` = CIRI dual arm (CRISPRa AND CRISPRi) |
| `threshold_a` | `-1` | CRISPRa UMI threshold; `-1` = auto KDE |
| `threshold_i` | `-1` | CRISPRi UMI threshold; `-1` = auto KDE |

### `ciri_step02_filter.R` (standalone)

| Argument | Default | Description |
|---|---|---|
| `--mito_lo/hi` | `0` / `10` | percent.mt bounds |
| `--ribo_lo` | `1` | percent.ribo lower bound |
| `--nGene_lo/hi` | `300` / `7000` | nFeature_RNA bounds |
| `--nUMI_lo` | `100` | nCount_RNA lower bound |
| `--num_dim` | `30` | PCA dimensions for Monocle3 |
| `--chosen_cds` | `cds_3` | Which of the 7 resolutions to finalise |

### `ciri_step06_trajectory()`

| Parameter | Default | Description |
|---|---|---|
| `clusters` | required | Cluster ID(s), e.g. `"5"` or `c("3","4")` |
| `root_gene` | required | Gene highest at trajectory start |
| `group` | required | Lineage name, e.g. `"muscle"` |
| `resolution` | `1e-3` | Sub-clustering resolution |
| `n_dims` | `50` | PCA dimensions for subset |

---

## Guides CSV format

No header. Three columns: `feature`, `type` (`a`/`i`), `fixed` (`f`/`v`).

```
MYOD_1,a,f          ← fixed CRISPRa (barcoding guide)
MYOD_2,a,f
NANOG,i,f           ← fixed CRISPRi (barcoding guide)
OCT4,i,f
SOX2,i,f
BAF60C_1,a,v        ← variable CRISPRa (perturbation guide)
BAF60C_2,a,v
CTCF_1,i,v          ← variable CRISPRi (perturbation guide)
NTCa,a,v            ← non-targeting control
NTCi,i,v
```

Use `guides_harmonised.csv` (from `harmonise_guide_names.R`) with Step 01,
not the original guides CSV — the names must match exactly what is in the H5.

---

## Key fixes vs original pipeline

| Issue | Fix |
|---|---|
| Mitochondrial genes survive protein-coding filter | `^MT-` genes explicitly removed after biomaRt filter |
| Guide name mismatch (H5 uses collapsed names) | `harmonise_guide_names.R` pre-step |
| `GetAssayData(slot=)` defunct in Seurat v5 | Updated to `layer=` throughout |
| `Remotes:` field forced monocle3 install | Removed; users install via `install_dependencies()` |
| NAMESPACE multi-line `importFrom` broke install | One `importFrom` per line |
