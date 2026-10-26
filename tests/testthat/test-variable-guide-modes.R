## Variable-guide calling protocols: Mode 1 (one guide per arm) and
## Mode 2 (A/B pair per arm, ported from CIRI_unified.R lines 388-470).

mk_long <- function(rows, type = "a") {
  data.frame(cell_barcode = rep(names(rows), lengths(rows)),
             feature      = unlist(lapply(rows, names), use.names = FALSE),
             umi          = unlist(rows, use.names = FALSE),
             type         = type,
             fixed        = "v",
             stringsAsFactors = FALSE)
}

long <- mk_long(list(
  # matching A/B pair, clear separation from the third guide
  cell1 = c(NTCa_1A = 20L, NTCa_1B = 18L, PAX7_1A = 1L),
  # top1 and top2 come from different genes -> rejected by the pairing filter
  cell2 = c(PAX7_1A = 20L, MYF6_1B = 18L, NTCa_1A = 1L),
  # correct pair but (top1 + top2) / top3 = 5, below min_ratio_pair
  cell3 = c(ID3_2A = 3L, ID3_2B = 2L, NTCa_1A = 1L),
  # exactly two variable guides: top3 undefined
  cell4 = c(MYF5_1A = 30L, MYF5_1B = 25L),
  # matching pair, passes every threshold
  cell5 = c(SMARCD3_1A = 10L, SMARCD3_1B = 9L, MYF6_1A = 1L)
))

test_that("Mode 2 calls matching A/B pairs and collapses them with a semicolon", {
  res <- CIRI:::.assign_var_guides_pair(long, fixed_genes = "MYOD", type_filter = "a")
  expect_setequal(res$cell_barcode, c("cell1", "cell5"))
  expect_equal(res$top_guide[res$cell_barcode == "cell1"], "NTCa_1A;NTCa_1B")
  expect_equal(res$top_guide[res$cell_barcode == "cell5"], "SMARCD3_1A;SMARCD3_1B")
})

test_that("Mode 2 rejects mismatched pairs and pairs below the ratio threshold", {
  res <- CIRI:::.assign_var_guides_pair(long, fixed_genes = "MYOD", type_filter = "a")
  expect_false("cell2" %in% res$cell_barcode)   # different genes
  expect_false("cell3" %in% res$cell_barcode)   # ratio 5 < 10
})

test_that("missing_top3 switches the fate of cells with exactly two guides", {
  drop   <- CIRI:::.assign_var_guides_pair(long, "MYOD", "a", missing_top3 = "drop")
  accept <- CIRI:::.assign_var_guides_pair(long, "MYOD", "a", missing_top3 = "accept")
  expect_false("cell4" %in% drop$cell_barcode)
  expect_true("cell4" %in% accept$cell_barcode)
  expect_equal(accept$top_guide[accept$cell_barcode == "cell4"], "MYF5_1A;MYF5_1B")
})

test_that("Mode 1 rejects an A/B-pair library, which is why Mode 2 exists", {
  # with two guides of the same gene at comparable UMI counts the Mode 1 ratio
  # top1/top2 is close to 1 and never reaches min_ratio = 5
  res <- CIRI:::.assign_var_guides(long, fixed_genes = "MYOD", type_filter = "a")
  expect_equal(nrow(res), 0L)
})

test_that("Mode 1 returns a single guide name on a one-guide library", {
  single <- mk_long(list(
    cell1 = c(NTCa_1A = 40L, PAX7_1A = 2L),
    cell2 = c(MYF6_1A = 30L, ID3_2A = 1L)
  ))
  res <- CIRI:::.assign_var_guides(single, fixed_genes = "MYOD", type_filter = "a")
  expect_setequal(res$cell_barcode, c("cell1", "cell2"))
  expect_true(all(!grepl(";", res$top_guide)))
  expect_equal(res$top_guide[res$cell_barcode == "cell1"], "NTCa_1A")
})

test_that("Mode 2 stops when the A/B replicates have been collapsed", {
  collapsed <- mk_long(list(
    cell1 = c(NTCa_1 = 20L, PAX7_1 = 1L),
    cell2 = c(MYF6_1 = 15L, ID3_2 = 2L)
  ))
  expect_error(
    CIRI:::.assign_var_guides_pair(collapsed, fixed_genes = "MYOD", type_filter = "a"),
    "requires guide names ending in A and B"
  )
})

test_that("both modes keep the cell_barcode + top_guide contract", {
  m1 <- CIRI:::.assign_var_guides(long, "MYOD", "a")
  m2 <- CIRI:::.assign_var_guides_pair(long, "MYOD", "a")
  expect_equal(names(m1), c("cell_barcode", "top_guide"))
  expect_equal(names(m2), c("cell_barcode", "top_guide"))
})
