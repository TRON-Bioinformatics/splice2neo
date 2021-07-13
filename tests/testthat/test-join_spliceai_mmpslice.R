test_that("join_spliceai_mmpslice works", {

  spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")
  df_raw <- parse_spliceai(spliceai_file)
  df <- format_spliceai(df_raw)
  annot_df <- annotate_spliceai_junction(df, toy_transcripts, toy_transcripts_gr)

  mmsplice_file <- system.file("extdata", "mmsplice_pred.csv", package = "splice2neo")
  mmsplice_df <- parse_mmsplice(mmsplice_file)
  mmsplice_df$pathogenicity <- "xx"
  mmsplice_df$effect <- "xx"
  mmsplice_df$efficiency <- "xx"
  junc_annot <- annotate_mmsplice(mmsplice_df, toy_transcripts)

  df_comb <- combine_mut_junc(annot_df, junc_annot)

  expect_true(nrow(df_comb) == 160)

})
