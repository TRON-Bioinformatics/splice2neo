test_that("parse_spliceai works on SpliceAI example file", {

  spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")
  df <- parse_spliceai(spliceai_file)

  expect_true(nrow(df) >= 10)

})

test_that("parse_spliceai works on a single-variant SpliceAI example file", {
  spliceai_file <- system.file(
    "extdata",
    "spliceai_output.one_variant.vcf",
    package = "splice2neo"
  )

  df <- parse_spliceai(spliceai_file)

  expect_true(nrow(df) == 1)
})

test_that("parse_spliceai works on a SpliceAI example file without variants", {
  spliceai_file <- system.file(
    "extdata",
    "spliceai_output.no_variant.vcf",
    package = "splice2neo"
  )

  df <- parse_spliceai(spliceai_file)

  expect_true(nrow(df) == 0)
})

test_that("parse_spliceai_thresh works on SpliceAI example file", {

  spliceai_file <- system.file("extdata", "spliceai_thresh_output.vcf", package = "splice2neo")
  df <- parse_spliceai_thresh(spliceai_file)

  expect_true(nrow(df) >= 10)

})

test_that(
  "parse_spliceai_thresh works on a single-variant SpliceAI example file", {
  spliceai_file <- system.file(
    "extdata",
    "spliceai_thresh_output.one_variant.vcf",
    package = "splice2neo"
  )

  df <- parse_spliceai_thresh(spliceai_file)

  expect_true(nrow(df) == 1)
})

test_that(
  "parse_spliceai_thresh works on a SpliceAI example file without variants", {
  spliceai_file <- system.file(
    "extdata",
    "spliceai_thresh_output.no_variant.vcf",
    package = "splice2neo"
  )

  df <- parse_spliceai_thresh(spliceai_file)

  expect_true(nrow(df) == 0)
})
