test_that("junc_to_gr works on toy data", {

  gr1 <- junc_to_gr("chr5:10-50:+")

  junc_id <- c("chr5:10-50:+", "chr_special1:5-20:-")
  gr2 <- junc_to_gr(junc_id)

  expect_equal(length(gr1), 1)
  expect_equal(length(gr2), 2)
  expect_equal(as.character(GenomeInfoDb::seqnames(gr2)), c("chr5", "chr_special1"))
})
