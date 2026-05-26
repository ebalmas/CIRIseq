# CIRI — Single-Cell CRISPRa/i Screen Analysis Pipeline

An R package for end-to-end analysis of CIRI screens: perturbation assignment,
QC filtering, Monocle3 preprocessing, cluster enrichment, pseudotime
statistics, and gene signature scoring.

---

## Installation

CIRI installs quickly because heavy bioinformatics packages (`Seurat`,
`monocle3`, `hdf5r`, `biomaRt`) are **not installed automatically**.
They are optional dependencies — each step checks for what it needs at
runtime and gives you a clear error with the install command if something
is missing.

### Step 1 — Install CIRI (fast, requires R ≥ 4.1)

```r
install.packages("devtools")
remotes::install_github("ebalmas/CIRIseq", ref = "V3")
```

That's it. Only lightweight CRAN packages (`dplyr`, `ggplot2`, etc.) are
installed automatically.

### Step 2 — Install the bioinformatics dependencies when you're ready

```r
library(CIRI)

# See what's installed and what's missing
check_dependencies()

# Install everything at once
install_dependencies()

# Or install only what you need for specific steps
install_dependencies(steps = "00")       # biomaRt (Ensembl reference)
install_dependencies(steps = "01")       # hdf5r (reading .h5 files)
install_dependencies(steps = "02")       # hdf5r + Seurat + data.table
install_dependencies(steps = "03-08")    # monocle3 + igraph
```

`check_dependencies()` shows exactly what's ready:

```
CIRI dependency status:
--------------------------------------------------
  [OK]    Step 00 — Ensembl reference         ready
  [MISS]  Step 01 — Guide assignment          missing: hdf5r
  [MISS]  Step 02 — Filter                    missing: hdf5r, Seurat
  [MISS]  Steps 03-08 — Monocle3              missing: monocle3, igraph
--------------------------------------------------
  1/4 step groups ready. Run install_dependencies() to install missing packages.
```

If you call a step before its dependency is installed, you get:

```
Error: Package 'monocle3' is required for this step but is not installed.
Install it with: BiocManager::install('monocle3')
```

### hdf5r system requirement

`hdf5r` requires the HDF5 C library on your system **before** `install_dependencies()`:

```bash
# macOS
brew install hdf5

# Ubuntu / Debian
sudo apt-get install libhdf5-dev

# Fedora / RHEL
sudo dnf install hdf5-devel
```

---

## Workflow overview

```
Output/                          ← dated run folders (never edit manually)
  <YYMMDD>_step01_assignment_AB011/
    csv/   plots/   stats/   R_objects/
    to_scratch/                  ← files ready for the next step

scratch/                         ← staging area YOU control
  annotation_data.csv            ← only here after ciri_promote_scratch()
```

**The key idea — deliberate checkpoints:**

```r
# 1. Run a step
ciri_step01_assignment(data_dir = "/data/AB011", matrix = "matrix.h5", sample = "AB011")

# 2. Inspect the output in Output/<date>_step01_assignment_AB011/
#    check plots/, stats/, csv/

# 3. When happy, promote to scratch/ so the next step can read it
ciri_promote_scratch("step01_assignment", sample = "AB011")

# 4. Run the next step — reads automatically from scratch/
ciri_step02_filter(data_dir = "/data/AB011", matrix = "matrix.h5", sample = "AB011")
```

Nothing moves to `scratch/` automatically. You decide when you are happy.

---

## Quick start

```r
library(CIRI)

# Copy the analysis template to your working directory
file.copy(
  system.file("run_analysis_template.R", package = "CIRI"),
  "run_analysis_AB011.R"
)
```

Open `run_analysis_AB011.R`, fill in your paths, and run each block one at a time.

If you are working directly from the cloned repo (no install), source locally instead:

```r
invisible(lapply(list.files("R", pattern = "\\.R$", full.names = TRUE), source))
```

---

## Functions

### Pipeline steps

| Function | Reads from | Writes to |
|---|---|---|
| `ciri_step00_download_ref()` | Ensembl API | `to_scratch/ensembl_protein_coding_genes.csv` |
| `ciri_step01_assignment()` | `data_dir/` (H5 + guides.csv) | `to_scratch/annotation_data.csv` |
| `ciri_step02_filter()` | `scratch/` + `data_dir/` (H5) | `to_scratch/annotated_matrix.csv` |
| `ciri_step03_load()` | `scratch/` | `to_scratch/processed_cds.RData` |
| `ciri_step04_validation()` | `scratch/` | plots + csv only (terminal) |
| `ciri_step05_enrichment()` | `scratch/` | plots + csv only (terminal) |
| `ciri_step06_trajectory()` | `scratch/` | `to_scratch/processed_cds_<group>.RData` + `pseudotime_<group>.csv` |
| `ciri_step07_pseudotime()` | `scratch/` | plots + csv only (terminal) |
| `ciri_step08_signatures()` | `scratch/` | plots + csv only (terminal) |

### Scratch management

| Function | Description |
|---|---|
| `ciri_promote_scratch()` | Copy a step's `to_scratch/` → `scratch/` (your checkpoint) |
| `ciri_list_scratch()` | List files ready in a step's `to_scratch/` |
| `ciri_scratch_status()` | Show current `scratch/` contents |
| `ciri_clear_scratch()` | Clear `scratch/` before a fresh run (requires `confirm = TRUE`) |

---

## Output structure

Each step writes to a dated folder:

```
Output/<YYMMDD>_<step>_<sample>/
  csv/          tables (.csv)
  plots/        figures (.pdf)
  stats/        text summaries / run logs
  R_objects/    .RData objects (permanent archive)
  to_scratch/   files consumed by the next step
```

`to_scratch/` is populated automatically. Use `ciri_promote_scratch()` to move
files into `scratch/` when you are satisfied with the output. Nothing in
`scratch/` is changed without your explicit action.

Add `Output/` and `scratch/` to `.gitignore`. Commit `R/`, `inst/`,
`DESCRIPTION`, `NAMESPACE`, and your filled `run_analysis_<experiment>.R` files.

---

## Parameter reference

### `ciri_step01_assignment()`
| Parameter | Default | Description |
|---|---|---|
| `data_dir` | required | Folder with H5 matrix and guides.csv |
| `matrix` | required | H5 filename |
| `sample` | `"CIRI"` | Experiment name (used in output folder) |
| `output_root` | `"Output"` | Top-level output directory |
| `guides` | `"guides.csv"` | Guide annotation CSV (no header: feature, type, fixed) |
| `strategy` | `1` | `1` = single variable guide; `2` = dual variable guides |
| `threshold_a` | `-1` | CRISPRa UMI threshold; `-1` = auto-detect via KDE valley |
| `threshold_i` | `-1` | CRISPRi UMI threshold; `-1` = auto-detect via KDE valley |

### `ciri_step02_filter()`
| Parameter | Default | Description |
|---|---|---|
| `data_dir` | required | Folder with H5 matrix |
| `matrix` | required | H5 filename |
| `sample` | `"CIRI"` | Experiment name |
| `scratch_dir` | `"scratch"` | Staging area (populated by `ciri_promote_scratch()`) |
| `gene_ref` | `NULL` | Path to `ensembl_protein_coding_genes.csv`; if `NULL`, resolved from `scratch/` |
| `min_genes` | `250` | Min detected genes per cell |
| `min_umis` | `3` | Min total UMIs per gene |
| `remove_mt` | `TRUE` | Remove mitochondrial genes (`^MT-`) |
| `remove_rb` | `TRUE` | Remove ribosomal genes (`^RPS\|^RPL`) |

### `ciri_step03_load()`
| Parameter | Default | Description |
|---|---|---|
| `sample` | `"CIRI"` | Experiment name |
| `scratch_dir` | `"scratch"` | Staging area |
| `resolution` | `5e-5` | Leiden clustering resolution (increase → more clusters) |
| `n_dims` | `100` | PCA dimensions |
| `seed` | `1234597698` | Random seed |

### `ciri_step06_trajectory()`
| Parameter | Default | Description |
|---|---|---|
| `clusters` | required | Cluster ID(s) to subset, character vector, e.g. `"5"` or `c("3","4")` |
| `root_gene` | required | Gene marking the trajectory start (highest-expression node) |
| `group` | required | Short lineage name used in output filenames, e.g. `"muscle"` |
| `resolution` | `1e-3` | Re-clustering resolution within the subset |
| `n_dims` | `50` | PCA dimensions for the subset |
| `seed` | `42` | Random seed |

### `ciri_step07_pseudotime()`
| Parameter | Default | Description |
|---|---|---|
| `group` | required | Group name matching `ciri_step06_trajectory()` |
| `control` | required | Control `gene_comb` string, e.g. `"NTCa-NA"` |
| `min_cells` | `8` | Min cells per perturbation to run KS test |
| `run_per_sample` | `FALSE` | Also run analysis separately per sample |
| `ecdf_top_n` | `10` | Number of top hits to save as ECDF plots |

### `ciri_step08_signatures()`
| Parameter | Default | Description |
|---|---|---|
| `group` | required | Group name matching `ciri_step06_trajectory()` |
| `signatures` | `NULL` | Named list of gene vectors; `NULL` uses built-in myogenic signatures |

---

## Guides CSV format

No header. Three columns: `feature`, `type` (`a` or `i`), `fixed` (`f` = fixed guide, `v` = variable guide).

```
NTCa_1A,a,f
NTCa_1B,a,f
NTCi_1A,i,f
NTCi_1B,i,f
SOX2_g1,a,v
SOX2_g2,a,v
CTCF_g1,i,v
CTCF_g2,i,v
```

---

## Key fixes vs original pipeline

**Mitochondrial gene removal:** MT-encoded genes (e.g. `MT-CO1`) are classified as
`protein_coding` by Ensembl biotype and survive the `biomaRt` filter in the original
`anno_filter.R`. `ciri_step02_filter()` explicitly removes `^MT-` genes *after*
the protein-coding filter, controlled by `remove_mt = TRUE/FALSE`.

**Duplicate protein-coding filter call** in the original `anno_filter.R` collapsed to one.

**Ensembl connectivity:** `ciri_step00_download_ref()` tries four mirrors in sequence
(`www`, `useast`, `uswest`, `asia`) instead of one, and also downloads
`chromosome_name` to flag mitochondrial genes in the reference.
