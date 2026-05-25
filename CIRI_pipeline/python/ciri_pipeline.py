#!/usr/bin/env python3
"""
CIRI Pipeline Runner
====================
Runs the CIRI single-cell CRISPR screen analysis pipeline locally
(no Docker required).

Steps:
  00  Download Ensembl protein-coding gene reference (run once)
  01  Perturbation assignment
  02  Annotation & filtering
  03  Data loading & preprocessing

Usage
-----
  # Run full pipeline:
  python ciri_pipeline.py run \
      --dir       /path/to/data \
      --matrix    filtered_feature_bc_matrix.h5 \
      --strategy  1

  # Run a specific step:
  python ciri_pipeline.py step 02 \
      --dir    /path/to/data \
      --matrix filtered_feature_bc_matrix.h5

  # Download gene reference (once):
  python ciri_pipeline.py download-ref --out ensembl_protein_coding_genes.csv

  # Check R dependencies:
  python ciri_pipeline.py check-deps
"""

import argparse
import logging
import os
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

logging.basicConfig(
    level=logging.INFO,
    format="[%(levelname)s] %(message)s"
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Locate the R scripts directory (same dir as this script)
# ---------------------------------------------------------------------------
SCRIPTS_DIR = Path(__file__).resolve().parent.parent / "R"


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------
@dataclass
class PipelineConfig:
    data_dir:        str
    matrix:          str
    guides:          str   = "guides.csv"
    strategy:        int   = 1
    threshold_a:     float = -1.0
    threshold_i:     float = -1.0
    gene_ref:        str   = "ensembl_protein_coding_genes.csv"
    min_genes:       int   = 250
    min_umis:        int   = 3
    remove_mt:       bool  = True
    remove_rb:       bool  = True
    resolution:      float = 5e-5
    seed:            int   = 1234597698
    rscript:         str   = "Rscript"

    def as_dict(self):
        return {k: v for k, v in self.__dict__.items()}


@dataclass
class StepResult:
    step:    str
    success: bool
    message: str = ""


# ---------------------------------------------------------------------------
# R script runner
# ---------------------------------------------------------------------------
def run_rscript(script: Path, args: List[str], rscript: str = "Rscript") -> StepResult:
    """Run an R script and return a StepResult."""
    cmd = [rscript, str(script)] + args
    log.info("Running: %s", " ".join(cmd))

    result = subprocess.run(cmd, capture_output=False, text=True)

    if result.returncode != 0:
        return StepResult(
            step    = script.name,
            success = False,
            message = f"R script exited with code {result.returncode}"
        )
    return StepResult(step=script.name, success=True)


# ---------------------------------------------------------------------------
# Individual steps
# ---------------------------------------------------------------------------
def step_download_ref(out_path: str, rscript: str = "Rscript") -> StepResult:
    script = SCRIPTS_DIR / "00_download_gene_ref.R"
    return run_rscript(script, [out_path], rscript=rscript)


def step_perturbation_assignment(cfg: PipelineConfig) -> StepResult:
    script = SCRIPTS_DIR / "01_perturbation_assignment.R"
    args = [
        "--dir",         cfg.data_dir,
        "--matrix",      cfg.matrix,
        "--guides",      cfg.guides,
        "--strategy",    str(cfg.strategy),
        "--threshold_a", str(cfg.threshold_a),
        "--threshold_i", str(cfg.threshold_i),
    ]
    return run_rscript(script, args, rscript=cfg.rscript)


def step_anno_filter(cfg: PipelineConfig) -> StepResult:
    script = SCRIPTS_DIR / "02_anno_filter.R"
    args = [
        "--dir",       cfg.data_dir,
        "--matrix",    cfg.matrix,
        "--gene_ref",  cfg.gene_ref,
        "--min_genes", str(cfg.min_genes),
        "--min_umis",  str(cfg.min_umis),
        "--remove_mt", str(cfg.remove_mt),
        "--remove_rb", str(cfg.remove_rb),
    ]
    return run_rscript(script, args, rscript=cfg.rscript)


def step_load(cfg: PipelineConfig) -> StepResult:
    script = SCRIPTS_DIR / "03_load.R"
    args = [
        "--dir",        cfg.data_dir,
        "--matrix",     "annotated_matrix.csv",
        "--resolution", str(cfg.resolution),
        "--seed",       str(cfg.seed),
    ]
    return run_rscript(script, args, rscript=cfg.rscript)


STEPS = {
    "00": ("Download gene reference",    lambda cfg: step_download_ref(cfg.gene_ref, cfg.rscript)),
    "01": ("Perturbation assignment",     step_perturbation_assignment),
    "02": ("Annotation & filtering",      step_anno_filter),
    "03": ("Load & preprocess",           step_load),
}


# ---------------------------------------------------------------------------
# Dependency checker
# ---------------------------------------------------------------------------
R_PACKAGES = [
    "dplyr", "tidyr", "stringr", "purrr", "ggplot2",
    "quantmod", "pracma", "hdf5r", "Matrix", "zoo",
    "scales", "Seurat", "biomaRt", "data.table",
    "monocle3", "viridis", "ggrepel", "gtools",
]

def check_deps(rscript: str = "Rscript") -> bool:
    """Check that all required R packages are installed."""
    pkgs_str = ", ".join(f'"{p}"' for p in R_PACKAGES)
    code = f"""
pkgs <- c({pkgs_str})
missing <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(missing) > 0) {{
  cat("MISSING:", paste(missing, collapse=", "), "\\n")
  quit(status = 1)
}} else {{
  cat("All packages available.\\n")
}}
"""
    result = subprocess.run(
        [rscript, "-e", code],
        capture_output=True, text=True
    )
    print(result.stdout.strip())
    if result.returncode != 0:
        log.error("Missing packages: %s", result.stdout.strip())
        return False
    return True


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="CIRI Analysis Pipeline (local, no Docker)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ---- run (all steps) ---------------------------------------------------
    run_p = sub.add_parser("run", help="Run the full pipeline")
    _add_pipeline_args(run_p)

    # ---- step (single step) ------------------------------------------------
    step_p = sub.add_parser("step", help="Run a single step")
    step_p.add_argument("step_id", choices=list(STEPS.keys()),
                        help="Step number (00, 01, 02, 03)")
    _add_pipeline_args(step_p)

    # ---- download-ref ------------------------------------------------------
    ref_p = sub.add_parser("download-ref", help="Download Ensembl gene reference")
    ref_p.add_argument("--out",     default="ensembl_protein_coding_genes.csv")
    ref_p.add_argument("--rscript", default="Rscript")

    # ---- check-deps --------------------------------------------------------
    dep_p = sub.add_parser("check-deps", help="Check R package dependencies")
    dep_p.add_argument("--rscript", default="Rscript")

    return parser


def _add_pipeline_args(p: argparse.ArgumentParser):
    # Required
    p.add_argument("--dir",         required=True,  help="Data directory")
    p.add_argument("--matrix",      required=True,  help="H5 matrix filename")
    # Optional
    p.add_argument("--guides",      default="guides.csv")
    p.add_argument("--strategy",    type=int,   default=1,    choices=[1, 2])
    p.add_argument("--threshold_a", type=float, default=-1.0,
                   help="CRISPRa UMI threshold (-1 = auto)")
    p.add_argument("--threshold_i", type=float, default=-1.0,
                   help="CRISPRi UMI threshold (-1 = auto)")
    p.add_argument("--gene_ref",    default="ensembl_protein_coding_genes.csv",
                   help="Path to ensembl_protein_coding_genes.csv")
    p.add_argument("--min_genes",   type=int,   default=250)
    p.add_argument("--min_umis",    type=int,   default=3)
    p.add_argument("--remove_mt",   type=lambda x: x.lower() != "false", default=True,
                   help="Remove mitochondrial genes (default: True)")
    p.add_argument("--remove_rb",   type=lambda x: x.lower() != "false", default=True,
                   help="Remove ribosomal genes (default: True)")
    p.add_argument("--resolution",  type=float, default=5e-5,
                   help="Monocle3 clustering resolution (default: 5e-5)")
    p.add_argument("--seed",        type=int,   default=1234597698)
    p.add_argument("--rscript",     default="Rscript",
                   help="Path to Rscript binary (default: Rscript)")


def args_to_config(args: argparse.Namespace) -> PipelineConfig:
    return PipelineConfig(
        data_dir    = str(Path(args.dir).resolve()),
        matrix      = args.matrix,
        guides      = args.guides,
        strategy    = args.strategy,
        threshold_a = args.threshold_a,
        threshold_i = args.threshold_i,
        gene_ref    = args.gene_ref,
        min_genes   = args.min_genes,
        min_umis    = args.min_umis,
        remove_mt   = args.remove_mt,
        remove_rb   = args.remove_rb,
        resolution  = args.resolution,
        seed        = args.seed,
        rscript     = args.rscript,
    )


def run_pipeline(cfg: PipelineConfig, step_ids: Optional[List[str]] = None):
    ids = step_ids or list(STEPS.keys())
    # Skip step 00 unless explicitly requested (gene ref should already exist)
    if step_ids is None:
        ids = [k for k in ids if k != "00"]

    results = []
    for sid in ids:
        label, fn = STEPS[sid]
        log.info("=== Step %s: %s ===", sid, label)
        res = fn(cfg)
        results.append(res)
        if not res.success:
            log.error("Step %s failed: %s", sid, res.message)
            log.error("Pipeline aborted.")
            sys.exit(1)
        log.info("Step %s completed successfully.", sid)

    log.info("Pipeline finished. %d steps completed.", len(results))


def main():
    parser = build_parser()
    args   = parser.parse_args()

    if args.command == "check-deps":
        ok = check_deps(rscript=args.rscript)
        sys.exit(0 if ok else 1)

    if args.command == "download-ref":
        cfg = PipelineConfig(data_dir=".", matrix="")
        cfg.rscript  = args.rscript
        cfg.gene_ref = args.out
        res = step_download_ref(args.out, rscript=args.rscript)
        if not res.success:
            log.error("Reference download failed: %s", res.message)
            sys.exit(1)
        return

    cfg = args_to_config(args)

    if args.command == "run":
        run_pipeline(cfg)
    elif args.command == "step":
        run_pipeline(cfg, step_ids=[args.step_id])


if __name__ == "__main__":
    main()
