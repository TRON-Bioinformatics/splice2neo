test_that("format_pangolin works on pangolin example file", {

  pangolin_file <- system.file("extdata", "spliceai_output.pangolin.vcf", package = "splice2neo")
  pangolin_df <- parse_pangolin(pangolin_file)

  df <- format_pangolin(pangolin_df)

  expect_true(nrow(df) >= 10)
  expect_true(all(df$score > 0))

})
