test_that("annotate_spliceai_junction works on toy example", {

  spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")
  df_raw <- parse_spliceai(spliceai_file)
  df <- format_spliceai(df_raw)

  annot_df <- annotate_spliceai_junction(df, toy_transcripts, toy_transcripts_gr)

  expect_true(nrow(annot_df) >= nrow(df))

})
