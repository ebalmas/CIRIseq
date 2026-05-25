# CIRI — Single-Cell CRISPRa/i Screen Analysis Pipeline

An R package for end-to-end analysis of CIRI screens: perturbation assignment,
QC filtering, Monocle3 preprocessing, cluster enrichment, pseudotime
statistics, and gene signature scoring.

---

## Installation

```r
# Install dependencies
install.packages("devtools")
BiocManager::install(c("biomaRt", "SingleCellExperiment"))
remotes::install_github("cole-trapnell-lab/monocle3")

# Install CIRI
devtools::install_github("yourusername/CIRI")
```

---

## Workflow overview

```
Output/                          ← dated run folders (never edit manually)
  <YYMMDD>_step01_assignment_AB011/
    csv/        plots/  stats/  R_objects/
    to_scratch/ ← files ready for the next step

scratch/                         ← staging area YOU control
  annotation_data.csv            ← copied here by ciri_promote_scratch()
  ...
```

**The key idea — deliberate checkpoints:**

```r
# 1. Run a step
ciri_step01_assignment(data_dir = "/data/AB011", matrix = "matrix.h5", sample = "AB011")

# 2. Inspect the output in Output/<date>_step01.../
#    plots/, stats/, csv/ — check everything looks right

# 3. When happy, promote to scratch/ so the next step can read it
ciri_promote_scratch("step01_assignment", sample = "AB011")

# 4. Run the next step — it reads from scratch/
ciri_step02_filter(data_dir = "/data/AB011", matrix = "matrix.h5", sample = "AB011")
```

Nothing moves to `scratch/` automatically. You decide when you are happy.

---

## Quick start

```r
library(CIRI)

# Copy the template and fill in your paths
file.copy(
  system.file("run_analysis_template.R", package = "CIRI"),
  "run_analysis_AB011.R"
)
```

Then open `run_analysis_AB011.R` and run each block interactively.

---

## Functions

| Function | Description |
|---|---|
| `ciri_step00_download_ref()` | Download Ensembl protein-coding gene reference (once) |
| `ciri_step01_assignment()` | Assign CRISPR guide identities to cells |
| `ciri_step02_filter()` | QC filtering + annotate expression matrix |
| `ciri_step03_load()` | Monocle3 CDS, UMAP, clustering |
| `ciri_step04_validation()` | Target knockdown / activation validation |
| `ciri_step05_enrichment()` | Cluster enrichment analysis |
| `ciri_step06_trajectory()` | Subclustering + pseudotime trajectory |
| `ciri_step07_pseudotime()` | KS test comparing pseudotime distributions |
| `ciri_step08_signatures()` | Gene set signature scoring |
| `ciri_promote_scratch()` | Copy `to_scratch/` → `scratch/` (your checkpoint) |
| `ciri_list_scratch()` | List files ready in a step's `to_scratch/` |
| `ciri_scratch_status()` | Show current `scratch/` contents |
| `ciri_clear_scratch()` | Clear `scratch/` before a fresh run |

---

## Output structure

Each step writes to a dated folder:

```
Output/<YYMMDD>_<step>_<sample>/
  csv/          tables (.csv)
  plots/        figures (.pdf)
  stats/        text summaries / logs
  R_objects/    .RData objects (permanent archive)
  to_scratch/   files needed by the next step
```

`to_scratch/` is read-only from your perspective — use `ciri_promote_scratch()`
to move files into `scratch/` when you are satisfied with the output.

Add `Output/` and `scratch/` to `.gitignore`. Commit only `R/`, the template,
and your filled `run_analysis_<experiment>.R`.

---

## Key fixes vs original pipeline

**Mitochondrial gene removal:** MT-encoded genes (e.g. `MT-CO1`) are protein-coding
by Ensembl biotype and survive the `biomaRt` filter. `ciri_step02_filter()` now
explicitly removes `^MT-` genes *after* the protein-coding filter (controlled by
`remove_mt = TRUE/FALSE`).

**Duplicate protein-coding filter call** in the original `anno_filter.R` removed.

**Ensembl connectivity:** `ciri_step00_download_ref()` tries four mirrors in
sequence and downloads `chromosome_name` to flag mito genes for reference.
