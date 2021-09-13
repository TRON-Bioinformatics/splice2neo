test_that("junc_to_cts works on example data with tx_id", {


  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  junc_id <- toy_junc_id[c(1, 6, 10)]
  tx_id <- toy_junc_id_enst[c(1, 6, 10)]
  transcripts <- toy_transcripts
  size = 400

  cts_df <- junc_to_cts(junc_id, transcripts, tx_id = tx_id, size = size, bsg = bsg)

  expect_true(nrow(cts_df) >= length(junc_id))
  expect_true(all(c("cts_seq", "cts_junc_pos", "cts_id") %in% names(cts_df)))
})

test_that("junc_to_cts works on example data without tx_id as input", {


  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  junc_id <- toy_junc_id[c(1, 6, 10)]
  transcripts <- toy_transcripts
  size = 400

  cts_df <- junc_to_cts(junc_id, transcripts, size = size, bsg = bsg)

  expect_true(nrow(cts_df) >= length(junc_id))
  expect_true(all(c("cts_seq", "cts_junc_pos", "cts_id") %in% names(cts_df)))
})
