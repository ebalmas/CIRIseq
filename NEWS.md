# CIRI 0.3.0

## Variable-guide calling protocols (Mode 1 and Mode 2)

`ciri_step01_assignment()` gains `variable_guides`, which selects how many
variable guides per arm the library carries. It is independent of `strategy`,
which continues to control which arms a cell must be positive for. The two were
previously conflated under one name.

* `variable_guides = 1` (default) — Mode 1, one variable guide per arm. The top
  guide is called when `top1 >= min_umi` and `top1 / top2 >= min_ratio`.
  This is the behaviour of all previous versions.
* `variable_guides = 2` — Mode 2, two variable guides per arm, the A/B pair of
  the same gene. Ported from the "Dual Guide Logic" branch of `CIRI_unified.R`
  in `alessandro-bertero/CIRI_bioinfo` (lines 388-470): the top three guides are
  ranked, the pair is called when `top1 + top2 >= min_umi_pair` and
  `(top1 + top2) / top3 >= min_ratio_pair`, top1 and top2 must be the A and B
  guide of the same gene and number, and the two names are returned collapsed
  with a semicolon (`NTCa_1A;NTCa_1B`).

Verified against the reference implementation on a 2-guide dataset: identical
cell sets and identical calls in both arms (6752 CRISPRa and 16 043 CRISPRi
cells, zero disagreements).

Mode 2 requires an input in which the A and B guides are distinct features. It
stops with an explanatory error when the guides CSV has no A/B names (checked
before the matrix is loaded) or when the matrix itself has them collapsed.

New thresholds, all at the values of the original pipeline: `min_umi`,
`min_ratio` (Mode 1), `min_umi_pair`, `min_ratio_pair` (Mode 2), plus
`missing_top3` for cells with exactly two variable guides in an arm
(`"drop"` reproduces published results, `"accept"` keeps them).

## Fixes

* `ciri_step01_assignment()` now accepts a guides CSV with or without a header
  row. `ciri_harmonise_guides()` writes `guides_harmonised.csv` *with* a header
  while step 01 read it with `header = FALSE`, which turned the header line into
  a guide called "feature" and produced a spurious "guide not found in data"
  warning on every run of the documented workflow.
* The `strategy = 2` documentation described a behaviour the code does not have.
  `.assign_dual()` returns dual cells together with CRISPRa-only and
  CRISPRi-only cells; the roxygen now says so.

## Known behaviour, unchanged

With H5 input the CRISPR table is dense (one row per cell per guide, zeros
included), exactly as in the original pipeline, so `top3` is usually present and
often zero. In Mode 2 this makes `(top1 + top2) / top3` infinite and leaves only
`top1 + top2 >= min_umi_pair` in force, and lets a pair be called when the
partner guide has zero UMI. Filtering `umi > 0` changes the result materially;
the published figures correspond to the unfiltered H5 path.
