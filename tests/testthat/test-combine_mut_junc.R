test_that("combine_mut_junc works", {

  # get spliceAI junctions
  spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")

  spliceai_annot_df <- parse_spliceai(spliceai_file) %>%
    format_spliceai() %>%
    annotate_mut_effect(toy_transcripts, toy_transcripts_gr)

  # get pangolin junctions
  pangolin_file <- system.file("extdata", "spliceai_output.pangolin.vcf", package = "splice2neo")

  pangolin_annot_df <- parse_pangolin(pangolin_file) %>%
    format_pangolin() %>%
    annotate_mut_effect(toy_transcripts, toy_transcripts_gr)

  # get mmsplice junctions
  mmsplice_file <- system.file("extdata", "mmsplice_pred.csv", package = "splice2neo")
  mmsplice_df <- parse_mmsplice(mmsplice_file)
  mmsplice_annot_df <- mmsplice_df %>%
    annotate_mmsplice(toy_transcripts)

  mmsplice_annot_df$pathogenicity <- "xx"
  mmsplice_annot_df$effect <- "xx"
  mmsplice_annot_df$efficiency <- "xx"

  junc_data_list = list(
    "spliceai" = spliceai_annot_df,
    "pangolin" = pangolin_annot_df,
    "mmsplice" = mmsplice_annot_df
  )

  n_unique <- bind_rows(junc_data_list) %>% select(mut_id, junc_id, tx_id) %>% n_distinct()

  df_comb <- combine_mut_junc(junc_data_list)

  expect_true(nrow(df_comb) >= 1)

})

test_that("combine_mut_junc works if one of the input table is empty", {

  # get spliceAI junctions
  spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")

  spliceai_annot_df <- parse_spliceai(spliceai_file) %>%
    format_spliceai() %>%
    annotate_mut_effect(toy_transcripts, toy_transcripts_gr)

  # get pangolin junctions
  pangolin_file <- system.file("extdata", "spliceai_output.pangolin.vcf", package = "splice2neo")

  pangolin_annot_df <- parse_pangolin(pangolin_file) %>%
    format_pangolin() %>%
    annotate_mut_effect(toy_transcripts, toy_transcripts_gr)

  # get mmsplice junctions
  mmsplice_file <- system.file("extdata", "mmsplice_pred.csv", package = "splice2neo")
  mmsplice_df <- parse_mmsplice(mmsplice_file)
  mmsplice_annot_df <- mmsplice_df %>%
    annotate_mmsplice(toy_transcripts)

  mmsplice_annot_df_empty <- mmsplice_annot_df %>%
    filter(row_number() < 1)

  junc_data_list = list(
    "spliceai" = spliceai_annot_df,
    "pangolin" = pangolin_annot_df,
    "mmsplice" = mmsplice_annot_df_empty
  )

  junc_data_list_wo_mmsplice = list(
    "spliceai" = spliceai_annot_df,
    "pangolin" = pangolin_annot_df
  )

  df_comb <- combine_mut_junc(junc_data_list)
  df_comb_wo_mmsplice <- combine_mut_junc(junc_data_list_wo_mmsplice)

  expect_true(nrow(df_comb) >= 1)
  expect_true(nrow(df_comb) == nrow(df_comb_wo_mmsplice))

  # column detected returned for all tools?
  expect_true(length(grep("detected", names(df_comb))) == length(junc_data_list))

})


test_that("combine_mut_junc provides expected number of rows", {

  # get spliceAI junctions
  spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")

  spliceai_annot_df <- parse_spliceai(spliceai_file) %>%
    format_spliceai() %>%
    annotate_mut_effect(toy_transcripts, toy_transcripts_gr)

  # get pangolin junctions
  pangolin_file <- system.file("extdata", "spliceai_output.pangolin.vcf", package = "splice2neo")

  pangolin_annot_df <- parse_pangolin(pangolin_file) %>%
    format_pangolin() %>%
    annotate_mut_effect(toy_transcripts, toy_transcripts_gr)


  merged_df <- spliceai_annot_df %>%
    dplyr::full_join(pangolin_annot_df, by = c("mut_id", "tx_id", "junc_id"), relationship = "many-to-many")

  junc_data_list_wo_mmsplice = list(
    "spliceai" = spliceai_annot_df,
    "pangolin" = pangolin_annot_df
  )

  df_comb <- combine_mut_junc(junc_data_list_wo_mmsplice)

  expect_true(nrow(df_comb) == nrow(merged_df))

})
