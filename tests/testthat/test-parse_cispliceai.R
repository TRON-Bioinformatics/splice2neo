test_that("parse_cispliceai_thresh works on CI-SpliceAI example file", {

  cispliceai_file <- system.file("extdata", "cispliceai_thresh_output.vcf", package = "splice2neo")
  df <- parse_cispliceai_thresh(cispliceai_file)

  expect_true(nrow(df) >= 10)

})

test_that("parse_cispliceai_thresh works on empty CI-SpliceAI example file", {

  cispliceai_file <- system.file("extdata", "cispliceai_thresh_output.empty.vcf", package = "splice2neo")
  df <- parse_cispliceai_thresh(cispliceai_file)

  expect_true(nrow(df) == 1)

})
