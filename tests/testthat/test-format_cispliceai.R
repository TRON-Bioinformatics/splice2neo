test_that("format_cispliceai_thresh works on CI-SpliceAI example file", {

  cispliceai_file <- system.file("extdata", "cispliceai_thresh_output.vcf", package = "splice2neo")
  df_raw <- parse_cispliceai_thresh(cispliceai_file)
  df <- format_cispliceai_thresh(df_raw)

  expect_true(nrow(df) >= 10)

  expected_columns <- c("mut_id", "effect", "score", "chr", "pos_rel", "pos")

  expect_equal( names(df), expected_columns )

})

test_that("format_cispliceai_thresh works on CI-SpliceAI example file with keep_gene_id option", {

  cispliceai_file <- system.file("extdata", "cispliceai_thresh_output.vcf", package = "splice2neo")
  df_raw <- parse_cispliceai_thresh(cispliceai_file)
  df <- format_cispliceai_thresh(df_raw)
  df_with_gene <- format_cispliceai_thresh(df_raw, transcripts_gr = toy_transcripts_gr )

  expect_equal(nrow(df), nrow(df_with_gene))

  expected_columns <- c("mut_id", "effect", "score", "chr", "pos_rel", "pos", "gene_id")

  expect_equal( names(df_with_gene), expected_columns )

})
