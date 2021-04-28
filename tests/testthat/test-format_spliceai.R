test_that("format_spliceai works on SpliceAI example file", {

  spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")
  df_raw <- parse_spliceai(spliceai_file)
  df <- format_spliceai(df_raw)

  expect_true(nrow(df) >= 10)

})
