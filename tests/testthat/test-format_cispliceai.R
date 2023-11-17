test_that("format_cispliceai_thresh works on CI-SpliceAI example file", {

  cispliceai_file <- system.file("extdata", "cispliceai_thresh_output.vcf", package = "splice2neo")
  df_raw <- parse_cispliceai_thresh(cispliceai_file)
  df <- format_cispliceai_thresh(df_raw)

  expect_true(nrow(df) >= 10)

})
