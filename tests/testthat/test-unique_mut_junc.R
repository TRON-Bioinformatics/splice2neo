test_that("unique_mut_junc works", {

  # get spliceAI junctions
  spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")
  splicai_df <- parse_spliceai(spliceai_file)

  spliceai_annot_df <- splicai_df %>%
    format_spliceai() %>%
    annotate_mut_effect(toy_transcripts, toy_transcripts_gr)

  spliceai_annot_df_unique <- unique_mut_junc(spliceai_annot_df)
  deduplicated_count <- spliceai_annot_df_unique %>% count(mut_id, junc_id, tx_id) %>% filter(n > 1)

  expect_true(nrow(spliceai_annot_df_unique) < nrow(spliceai_annot_df))
  expect_true(nrow(deduplicated_count) == 0)


})


test_that("unique_junc_mmsplice works", {

  mmsplice_annot_df <- tibble(
    mut_id = c("id1", "id1","id1", "id2", "id3"),
    junc_id = c("junc_id1", "junc_id1","junc_id1", "junc_id2", "junc_id3"),
    tx_id = c("tx_id1", "tx_id1","tx_id1", "tx_id2", "tx_id3"),
    event_type = c("exon inclusion", "exon inclusion", "exon inclusion","exon skipping", "exon inclusion"),
    delta_logit_psi = c(1, 2, 2,-1, 1),
  )


  mmsplice_annot_df_unique <- unique_junc_mmsplice(mmsplice_annot_df)
  deduplicated_count <- mmsplice_annot_df_unique %>% count(mut_id, junc_id, tx_id) %>% filter(n > 1)

  expect_true(nrow(mmsplice_annot_df_unique) < nrow(mmsplice_annot_df))
  expect_true(nrow(deduplicated_count) == 0)


})
