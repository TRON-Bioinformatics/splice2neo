test_that("add_context_seq works on toy example data", {

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  junc_df <- tibble::tibble(
    junc_id = toy_junc_id[c(1, 6, 10)]
  ) %>%
    add_tx(toy_transcripts)

  cts_df <- add_context_seq(junc_df, size = 400, bsg = bsg)

  expect_true(nrow(cts_df) == nrow(junc_df))
  expect_true(all(c("cts_seq", "cts_junc_pos", "cts_id") %in% names(cts_df)))
})
