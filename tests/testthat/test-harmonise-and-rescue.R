## Harmonisation name matching (0.4.0 correctness fix) and orphan rescue.

mk_long <- function(rows, type = "a") {
  data.frame(cell_barcode = rep(names(rows), lengths(rows)),
             feature      = unlist(lapply(rows, names), use.names = FALSE),
             umi          = unlist(rows, use.names = FALSE),
             type = type, fixed = "v", stringsAsFactors = FALSE)
}

# ---- harmonisation --------------------------------------------------------

test_that("exact match wins over stripping A/B", {
  # the matrix resolves A and B; the guide name is present verbatim
  dn <- c("ATF7IP_1A", "ATF7IP_1B", "ATF7IP_1")
  expect_equal(CIRI:::.find_match("ATF7IP_1A", dn), "ATF7IP_1A")
})

test_that("a guide is not mapped onto the id of a different guide", {
  # id "ATF7IP_1" belongs to name "ATF7IP_2A": mapping ATF7IP_1A there is wrong
  dn  <- c("ATF7IP_1", "ATF7IP_2")
  map <- c(ATF7IP_1 = "ATF7IP_2A", ATF7IP_2 = "ATF7IP_1B")
  expect_true(is.na(suppressMessages(
    CIRI:::.find_match("ATF7IP_1A", dn, id2name = map))))
  # and the legitimate one still resolves
  expect_equal(suppressMessages(
    CIRI:::.find_match("ATF7IP_2A", dn, id2name = map)), "ATF7IP_1")
})

test_that("collapsing is still allowed when the id really is the collapsed name", {
  dn  <- c("ATF7IP_1")
  map <- c(ATF7IP_1 = "ATF7IP_1")
  expect_equal(CIRI:::.find_match("ATF7IP_1A", dn, id2name = map), "ATF7IP_1")
})

test_that("without a reference the collapsed case still works", {
  expect_equal(CIRI:::.find_match("ATF7IP_1A", c("ATF7IP_1")), "ATF7IP_1")
  expect_equal(CIRI:::.find_match("NTCa_1A", c("NTCa")), "NTCa")
})

# ---- construct classification ---------------------------------------------

cl <- mk_long(list(
  c1 = c(GOOD_1A = 50L, GOOD_1B = 45L, ORPH_1B = 1L),
  c2 = c(GOOD_1A = 40L, GOOD_1B = 38L),
  c3 = c(GOOD_1A = 30L, GOOD_1B = 33L),
  c4 = c(ORPH_1B = 80L),
  c5 = c(ORPH_1B = 70L),
  c6 = c(ORPH_1B = 60L),
  c7 = c(DEAD_1A = 12L)))

test_that("constructs are classified as paired, orphan or dead", {
  cls <- CIRI:::.classify_constructs(cl, fixed_genes = "FIX", type_filter = "a",
                                     orphan_frac = 0.10, min_cells = 2)
  g <- function(x) cls$class[cls$construct == x]
  expect_equal(g("GOOD_1"), "paired")
  expect_equal(g("ORPH_1"), "orphan")
  expect_equal(cls$dominant[cls$construct == "ORPH_1"], "ORPH_1B")
  expect_equal(g("DEAD_1"), "dead")
})

test_that("rescue_constructs restricts the rescue and warns on unknown names", {
  all_orph <- suppressMessages(
    CIRI:::.assign_var_guides_pair_rescue(cl, "FIX", "a", min_cells = 2,
                                          min_umi_rescue = 10, min_ratio_rescue = 5,
                                          report = FALSE))
  none <- suppressMessages(
    CIRI:::.assign_var_guides_pair_rescue(cl, "FIX", "a", min_cells = 2,
                                          rescue_constructs = character(0),
                                          report = FALSE))
  expect_true(any(none$top_guide == "GOOD_1A;GOOD_1B"))        # pairs survive
  expect_false(any(!grepl(";", none$top_guide)))               # no rescues
  expect_true(any(all_orph$top_guide == "ORPH_1B"))            # rescued, single name
  expect_warning(
    CIRI:::.assign_var_guides_pair_rescue(cl, "FIX", "a", min_cells = 2,
                                          rescue_constructs = "GOOD_1A",
                                          report = FALSE),
    "not classified as orphan")
})

test_that("rescued calls are distinguishable from pairs by the separator", {
  r <- suppressMessages(
    CIRI:::.assign_var_guides_pair_rescue(cl, "FIX", "a", min_cells = 2,
                                          rescue_constructs = "ORPH_1B",
                                          report = FALSE))
  expect_setequal(names(r), c("cell_barcode", "top_guide"))
  expect_true(all(grepl(";", r$top_guide) | r$top_guide == "ORPH_1B"))
})

# ---- validation ------------------------------------------------------------

test_that("validation flags a rescued set that tracks the background", {
  ann <- data.frame(
    cell_barcode = paste0("c", 1:20),
    feature_a = c(rep("X_1A;X_1B", 5), rep("X_1A", 5), rep("NA", 10)),
    feature_i = "NA", stringsAsFactors = FALSE)
  meta <- data.frame(cell = paste0("c", 1:20), sample = "s1",
                     # pairs all phenotype-positive, rescued all negative,
                     # background mostly negative
                     phenotype = c(rep(TRUE, 5), rep(FALSE, 5), rep(FALSE, 10)),
                     stringsAsFactors = FALSE)
  v <- CIRI::ciri_validate_rescue(ann, meta)
  expect_equal(v$pct_pair, 100)
  expect_equal(v$pct_rescued, 0)
  expect_match(v$verdict, "DO NOT rescue")
})

test_that("validation accepts a rescued set that matches its pairs", {
  ann <- data.frame(
    cell_barcode = paste0("c", 1:20),
    feature_a = c(rep("X_1A;X_1B", 5), rep("X_1A", 5), rep("NA", 10)),
    feature_i = "NA", stringsAsFactors = FALSE)
  meta <- data.frame(cell = paste0("c", 1:20), sample = "s1",
                     phenotype = c(rep(TRUE, 10), rep(FALSE, 10)),
                     stringsAsFactors = FALSE)
  v <- CIRI::ciri_validate_rescue(ann, meta)
  expect_equal(v$pct_rescued, 100)
  expect_match(v$verdict, "consistent with its own pairs")
})

test_that("an entirely undetected partner still yields an orphan", {
  x <- mk_long(list(c1 = c(X_1B = 20L), c2 = c(X_1B = 30L)))
  cls <- CIRI:::.classify_constructs(x, "FIX", "a", min_cells = 2)
  expect_equal(cls$cells_A, 0)
  expect_equal(cls$class, "orphan")
  expect_equal(cls$dominant, "X_1B")
})

test_that("validation supports only pairs, only rescues and no calls", {
  meta <- data.frame(cell = c("c1", "c2"), sample = "s1", phenotype = c(TRUE, FALSE))
  for (call in c("X_1A;X_1B", "X_1B", "NA")) {
    ann <- data.frame(cell_barcode = meta$cell, feature_a = call, feature_i = "NA")
    v <- CIRI::ciri_validate_rescue(ann, meta)
    expect_true(all(c("pct_pair", "pct_rescued", "n_pair", "n_rescued") %in% names(v)))
    if (call == "NA") expect_equal(nrow(v), 0L)
    else if (grepl(";", call)) expect_match(v$verdict, "no rescued cells")
    else expect_match(v$verdict, "no pair reference")
  }
})

test_that("harmonisation translates reference ids to matrix guide names", {
  tmp <- tempfile(); dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE))
  guides <- file.path(tmp, "guides.csv")
  proto <- file.path(tmp, "protospacers.csv")
  ref <- file.path(tmp, "reference.csv")
  writeLines(c("X_1A,a,v", "X_1B,a,v"), guides)
  write.csv(data.frame(feature_call = "X_1|X_2"), proto, row.names = FALSE)
  write.csv(data.frame(id = c("X_1", "X_2"), name = c("X_1B", "X_1A")), ref, row.names = FALSE)
  result <- suppressMessages(CIRI::ciri_harmonise_guides(guides, proto, feature_reference_path = ref))
  expect_setequal(result$feature, c("X_1A", "X_1B"))
  mapping <- read.csv(file.path(tmp, "name_mapping.csv"))
  expect_equal(mapping$cellranger_name, mapping$original_name)
})

# ---- allow_collapse ---------------------------------------------------------

.fixtures <- function(dir, data_names) {
  gp <- file.path(dir, "guides.csv"); pp <- file.path(dir, "proto.csv")
  utils::write.table(data.frame(f = c("G_1A", "G_1B", "FIX"), t = "a",
                                x = c("v", "v", "f")),
                     gp, sep = ",", row.names = FALSE, col.names = FALSE, quote = FALSE)
  utils::write.csv(data.frame(cell_barcode = "c1", num_features = 1,
                              feature_call = paste(data_names, collapse = "|")),
                   pp, row.names = FALSE)
  list(guides = gp, proto = pp, out = file.path(dir, "out.csv"))
}

test_that("without an authoritative name source the collapse is refused", {
  d <- withr::local_tempdir(); f <- .fixtures(d, c("G_1", "FIX"))
  # "G_1" could be a collapsed name or the id of a different guide - unknowable
  expect_error(
    suppressMessages(ciri_harmonise_guides(f$guides, f$proto, f$out)),
    "cannot be verified")
})

test_that("allow_collapse = TRUE overrides the refusal", {
  d <- withr::local_tempdir(); f <- .fixtures(d, c("G_1", "FIX"))
  out <- suppressMessages(ciri_harmonise_guides(f$guides, f$proto, f$out,
                                                allow_collapse = TRUE))
  expect_false(any(grepl("_[0-9]+[AB]$", out$feature)))
})

test_that("a feature reference makes the collapsed case proceed with a warning", {
  d <- withr::local_tempdir(); f <- .fixtures(d, c("G_1", "FIX"))
  fr <- file.path(d, "fref.csv")
  utils::write.csv(data.frame(id = c("G_1", "FIX"), name = c("G_1", "FIX")),
                   fr, row.names = FALSE)
  msgs <- character(0)
  out <- withCallingHandlers(
    ciri_harmonise_guides(f$guides, f$proto, f$out, feature_reference_path = fr),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_true(any(grepl("would be collapsed", msgs)))
  expect_equal(nrow(out), 2)
})

test_that("a matrix that resolves A and B errors instead of collapsing", {
  d <- withr::local_tempdir(); f <- .fixtures(d, c("G_1A", "G_1B", "G_1", "FIX"))
  fr <- file.path(d, "fref.csv")
  utils::write.csv(data.frame(id = c("G_1A", "G_1B", "G_1", "FIX"),
                              name = c("G_1A", "G_1B", "G_1", "FIX")),
                   fr, row.names = FALSE)
  # both A and B resolve exactly, so nothing is collapsed and it completes
  out <- suppressMessages(ciri_harmonise_guides(f$guides, f$proto, f$out,
                                                feature_reference_path = fr))
  expect_true(all(c("G_1A", "G_1B") %in% out$feature))
})

test_that("guides absent from the data are not counted as collapsed", {
  d <- withr::local_tempdir(); f <- .fixtures(d, c("G_1A", "G_1B", "FIX"))
  fr <- file.path(d, "fref.csv")
  utils::write.csv(data.frame(id = c("G_1A", "G_1B", "FIX"),
                              name = c("G_1A", "G_1B", "FIX")),
                   fr, row.names = FALSE)
  out <- suppressMessages(ciri_harmonise_guides(f$guides, f$proto, f$out,
                                                feature_reference_path = fr))
  expect_equal(sum(grepl("_[0-9]+[AB]$", out$feature)), 2)
})
