test_that("is_canonical works on basic example", {

  ref_junc <- toy_junc_id[1:3]

  exon_gr <- unlist(toy_transcripts)

  can <- is_canonical(toy_junc_id, ref_junc, exon_gr)

  expect_true(all(can[1:3]))
})


