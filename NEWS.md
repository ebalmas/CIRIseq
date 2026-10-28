# CIRI 0.4.1

Two defects in the 0.4.0 collapse guard, both found by running the published
0.4.0 on AB014's real files rather than on fixtures.

* **The guard could not fire on the call that needed it.** `allow_collapse = NA`
  decided by looking for A/B names among the names being matched. Without
  `matrix_path` or `feature_reference_path` those names come from the
  protospacer CSV, i.e. the reference's `id` column, which for a multi-guide
  reference reads `GENE_1`, `GENE_2`, ... and contains no A/B at all. The guard
  concluded "already collapsed", warned, and collapsed anyway - the exact
  outcome 0.4.0 was written to prevent. On AB014 a bare
  `ciri_harmonise_guides(guides_2.csv, protospacer_calls_per_cell.csv)` still
  returned 30 rows with no A/B names left.

  `NA` now distinguishes three cases: A/B present among the matched names ->
  error; absent, with an authoritative source (`matrix_path` or
  `feature_reference_path`) -> warning and proceed, the genuinely collapsed
  reference; absent, with no authoritative source -> **error**, because the two
  cannot be told apart. `allow_collapse = TRUE` still overrides.

* **Undetected guides were counted as collapsed.** The guard compared the count
  of A/B names in the guide list against those surviving in the output, so a
  guide that matched *nothing* - a construct not detected in the experiment -
  counted as collapsed. On AB014 with `feature_reference_path` this produced a
  spurious error about "3 of 52" that referred to the three undetected MEF2C
  guides. The count is now taken over guides that matched onto a non-A/B name.

Verified on AB014: the bare call errors; `matrix_path` completes with 57 rows
and all 52 A/B names kept; `feature_reference_path` completes with 54 rows and
49 A/B names, the three missing being exactly the undetected MEF2C guides; a
genuinely collapsed reference still collapses with a warning.

# CIRI 0.4.0

## Guide-name harmonisation: correctness fix

`ciri_harmonise_guides()` matched guide names against the names found in
`crispr_analysis/protospacer_calls_per_cell.csv`. That file reports the **`id`**
column of the feature reference, while `filtered_feature_bc_matrix.h5` reports
the **`name`** column. When the two differ - a 4-guides-per-gene reference
numbers its rows `GENE_1` ... `GENE_4` under `id` while `name` carries the A/B
label - the rule that strips a trailing A/B always found a match, and mapped a
guide onto a different guide's feature.

Measured on AB014: `id` and `name` differ for 52 of 57 rows; of the 50 A/B
guides mapped onto an existing id, **40 landed on an id belonging to a different
guide**. The surviving rows were not merely collapsed, their guide-to-feature
correspondence was wrong.

* Exact match is now tried **first**. A guide name present in the matrix is
  never a candidate for collapsing.
* New `feature_reference_path` builds the `id`/`name` map and enables a guard
  that refuses to map a guide onto the id of a different guide.
* New `matrix_path` matches against `features/name` from the H5 itself, which is
  what the rest of the pipeline consumes. Takes precedence; needs hdf5r.
* New `allow_collapse`. By default harmonisation now **errors** when it would
  collapse A/B names on a matrix that resolves A and B separately, rather than
  silently producing a guide list on which `variable_guides = 2` cannot run.

Harmonisation is still required for the genuinely collapsed case - a reference
that declares one feature per construct. It no longer guesses.

## Rescue of constructs with an undetectable partner guide

The Mode 2 pairing filter confirms vector identity by requiring both guides of a
construct in the same cell. When one guide is undetectable across the whole
experiment, that confirmation is unavailable for every cell carrying it, and
requiring it discards real cells.

`ciri_step01_assignment()` gains `rescue_orphans`, `rescue_constructs`,
`orphan_frac`, `min_cells_construct`, `min_umi_rescue` and `min_ratio_rescue`.
Rescued calls carry ONE guide name where pairs carry two joined by `;`, so
`grepl(";", feature_a)` separates them and gene-level columns are unchanged.
`assignment_summary.txt` records the setting and the thresholds.

**Eligibility is not justification.** On AB014 four constructs were orphans and
only one survived validation: the rescued cells of `NTCa_1`, `JMJD3_2` and
`YBX1_1` tracked the sample background instead of the construct's own
pair-confirmed cells (NTCa pairs 95.7% / 83.3% mature vs rescued 10.5% / 3.2%),
i.e. contamination by an abundant barcode. Only `MEF2C_2B`, whose partner is
detected in zero cells, was enriched over background in both samples (67.6% vs
35.9%, 20.1% vs 4.1%). `rescue_constructs` should therefore always be given
explicitly; leaving it `NULL` rescues every orphan and on AB014 would have added
1387 wrong CRISPRa calls.

## Package validation

* Generate and ship `man/` help pages and reproducible namespace imports.
* Expose the name-matching helper internally so its regression tests run.
* Handle constructs with no detections of one partner and validation tables
  containing only pairs, only rescues, or no calls.
* Avoid reporting unvalidated constructs without pair references as validated.

## Rescue QC

`ciri_rescue_report()`, `ciri_validate_rescue()`,
`ciri_plot_rescue_validation()` and `ciri_rescue_qc()`. The validation cannot
live in step 01 because it needs a per-cell phenotype that only exists after
clustering; `ciri_rescue_qc()` writes `rescue_validation.csv` and
`rescue_validation.pdf` and stops if a construct named in `rescue_constructs`
is inconsistent with its own pairs. The phenotype is an argument, not fixed to
any particular cluster set.

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
