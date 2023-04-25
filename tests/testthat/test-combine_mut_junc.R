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
  mmsplice_annot_df <- parse_mmsplice(mmsplice_file) %>%
    annotate_mmsplice(toy_transcripts)

  mmsplice_annot_df$pathogenicity <- "xx"
  mmsplice_annot_df$effect <- "xx"
  mmsplice_annot_df$efficiency <- "xx"
  # junc_annot <- annotate_mmsplice(mmsplice_df, toy_transcripts)

  junc_data_list = list(
    "spliceai" = spliceai_annot_df,
    "pangolin" = pangolin_annot_df,
    "mmsplice" = mmsplice_annot_df
  )

  df_comb <- combine_mut_junc(junc_data_list)

  expect_true(nrow(df_comb) >= 1)

})


test_that("combine_mut_junc works if altered transcript was predicted with several tools", {

  # get spliceAI junctions
  spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")

  spliceai_annot_df <- tibble::tibble(
    mut_id = c("mut_id1","mut_id1", "mut_id2"),
    junc_id = c("junc_id1","junc_id1","junc_id2"),
    tx_id = c("tx_id1","tx_id1","tx_id2"),
    event_type = c("IR", "A3SS", "ES")
  )

  pangolin_annot_df <- tibble::tibble(
    mut_id = rep("mut_id1",2),
    junc_id = rep("junc_id1",2),
    tx_id = rep("tx_id1",2),
    event_type = c("IR", "A3SS")
  )

  mmsplice_annot_df <- tibble::tibble(
    mut_id = "mut_id2",
    junc_id = "junc_id2",
    tx_id = "tx_id2",
    event_type = "ES"
  )

  junc_data_list = list(
    "spliceai" = spliceai_annot_df,
    "pangolin" = pangolin_annot_df,
    "mmsplice" = mmsplice_annot_df
  )

  df_comb <- combine_mut_junc(junc_data_list)

  expect_true(nrow(df_comb) ==  3)

})

