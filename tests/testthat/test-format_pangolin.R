test_that("format_pangolin works on pangolin example file", {

  pangolin_file <- system.file("extdata", "spliceai_output.pangolin.vcf", package = "splice2neo")
  pangolin_df <- parse_pangolin(pangolin_file)

  df <- format_pangolin(pangolin_df)

  expect_true(nrow(df) >= 10)
  expect_true(!"gene_id" %in% names(df))
  expect_true(all(df$score > 0))

})

test_that("format_pangolin works on pangolin example file with keep_gene_id = TRUE", {

  pangolin_file <- system.file("extdata", "spliceai_output.pangolin.vcf", package = "splice2neo")
  pangolin_df <- parse_pangolin(pangolin_file)

  df <- format_pangolin(pangolin_df, keep_gene_id = TRUE)

  expect_true(nrow(df) >= 10)
  expect_true(all(df$score > 0))
  expect_true("gene_id" %in% names(df))
  expect_true(all(!is.na(df$gene_id)))

})
