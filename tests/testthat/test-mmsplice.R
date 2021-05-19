test_that("parse_mmsplice works with example file", {

  mmsplice_file <- system.file("extdata", "mmsplice_pred.csv", package = "splice2neo")
  df <- parse_mmsplice(mmsplice_file)

  expect_true(nrow(df) >= 10)

})


test_that("annotate_mmsplice() works with toy data", {

  mmsplice_file <- system.file("extdata", "mmsplice_pred.csv", package = "splice2neo")
  mmsplice_df <- parse_mmsplice(mmsplice_file)

  junc_annot <- annotate_mmsplice(mmsplice_df, toy_transcripts)

  expect_true(nrow(junc_annot) >= 10)


})
