test_that("junc_to_gr works on toy data", {

  gr1 <- junc_to_gr("chr1_5_10_+")

  junc_id <- c("chr1_5_10_+", "chr_special1_5_20_-")
  gr2 <- junc_to_gr(junc_id)

  expect_equal(length(gr1), 1)
  expect_equal(length(gr2), 2)
  expect_equal(as.character(GenomeInfoDb::seqnames(gr2)), c("chr1", "chr_special1"))
})
