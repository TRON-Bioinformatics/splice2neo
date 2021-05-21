test_that("junc_to_peptide works on example data", {

  # testthat::skip("Not implemented yet")

  toy_junc_id
  toy_cds

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  pep_df <- junc_to_peptide(toy_junc_id, toy_cds, size = 30, bsg = bsg)

  expect_true(nrow(pep_df) > 0)
  expect_true(nrow(pep_df) >= length(toy_junc_id))
  expect_equal(unique(pep_df$junc_id), unique(toy_junc_id))

  # check that peptides have expcted size
  expect_true(all(stringr::str_length(pep_df$peptide_context) <= 30, na.rm = TRUE))

})
