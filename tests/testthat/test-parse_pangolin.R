test_that("parse_pangolin works on provided example file", {

  pangolin_file <- system.file("extdata", "spliceai_output.pangolin.vcf", package = "splice2neo")
  pangolin_df <- parse_pangolin(pangolin_file)

  expect_true(nrow(pangolin_df) > 0)
  expect_true(ncol(pangolin_df) == 11)

  expect_true(all(
    c("pos_rel", "pangolin_score") %in% names(pangolin_df)
    ))

})

test_that("parse_pangolin works on one-variant example file", {

  pangolin_file <- system.file(
    "extdata",
    "spliceai_output.pangolin.one_variant.vcf",
    package = "splice2neo"
  )
  pangolin_df <- parse_pangolin(pangolin_file)

  # One variant, but two sites.
  expect_true(nrow(pangolin_df) == 2)
  expect_true(ncol(pangolin_df) == 11)

  expect_true(all(
    c("pos_rel", "pangolin_score") %in% names(pangolin_df)
    ))

})

test_that("parse_pangolin works on example file without variants", {

  pangolin_file <- system.file(
    "extdata",
    "spliceai_output.pangolin.no_variant.vcf",
    package = "splice2neo"
  )
  pangolin_df <- parse_pangolin(pangolin_file)

  expect_true(nrow(pangolin_df) == 0)
  expect_true(ncol(pangolin_df) == 11)

  expect_true(all(
    c("pos_rel", "pangolin_score") %in% names(pangolin_df)
    ))

})
