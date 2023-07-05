test_that("parse_pangolin works on provided example file", {

  pangolin_file <- system.file("extdata", "spliceai_output.pangolin.vcf", package = "splice2neo")
  pangolin_df <- parse_pangolin(pangolin_file)

  expect_true(nrow(pangolin_df) > 0)

  expect_true(all(
    c("pos_rel", "pangolin_score") %in% names(pangolin_df)
    ))

})
