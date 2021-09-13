test_that("junc_to_peptide works on example data", {

  # testthat::skip("Not implemented yet")

  # toy_junc_id
  # toy_cds

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  pep_df <- junc_to_peptide(toy_junc_id, toy_cds, size = 30, bsg = bsg)

  expect_true(nrow(pep_df) > 0)
  expect_true(nrow(pep_df) >= length(toy_junc_id))
  expect_equal(unique(pep_df$junc_id), unique(toy_junc_id))

  # check that peptides have expcted size
  expect_true(all(stringr::str_length(pep_df$peptide_context) <= 30, na.rm = TRUE))

})

test_that("junc_to_peptide works on example data with input tx ids", {

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  junc_id <- toy_junc_id[c(1, 6, 10)]
  tx_id <- toy_junc_id_enst[c(1, 6, 10)]
  cds <- toy_cds

  pep_df <- junc_to_peptide(junc_id, cds, tx_id = tx_id, size = 30, bsg = bsg)

  expect_true(nrow(pep_df) > 0)
  expect_true(nrow(pep_df) >= length(junc_id))
  expect_equal(unique(pep_df$junc_id), unique(junc_id))

  # check that peptides have expcted size
  expect_true(all(stringr::str_length(pep_df$peptide_context) <= 30, na.rm = TRUE))

})


test_that("seq_extract_nonstop works on example data with single seq", {

  seq <- "QIP*LGSNSLLFPYQLMAGSTRP*SWALGC"
  pos <- 14

  df <-  seq_extract_nonstop(seq, pos)

  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 1)
})

test_that("seq_extract_nonstop works on example data with multiple seq", {

  seq <- c(
   "QIP*LGSNSLLFPYQLMAGSTRP*SWALGC",
   "LKMRGDTNDILSHLD*REQRVGQ*AEAASP"
  )
  pos <- c(14, 14)
  df <- splice2neo:::seq_extract_nonstop(seq, pos)

  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 2)

})
