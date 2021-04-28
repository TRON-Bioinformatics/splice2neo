test_that("parse_spliceai works on SpliceAI example file", {

  spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")
  df <- parse_spliceai(spliceai_file)

  expect_true(nrow(df) >= 10)

})
