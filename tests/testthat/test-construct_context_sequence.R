test_that("map_junc_transcript works", {
  junc_id <- "chr2_100722033_100722195_+"
  junc_df <- tibble(junc_id = junc_id) %>%
    separate(
      junc_id,
      sep = "_",
      into = c("chr", "pos1", "pos2", "strand"),
      remove = FALSE
    ) %>%
    mutate(pos1 = as.integer(pos1), pos2 = as.integer(pos2))
  junc_pos1 <-
    GenomicRanges::GRanges(junc_df$chr, junc_df$pos1, strand = junc_df$strand)
  junc_pos2 <-
    GenomicRanges::GRanges(junc_df$chr, junc_df$pos2, strand = junc_df$strand)
  transcripts_covering_junction <-
    map_junc_transcript(
      junc_id = junc_id,
      junc_pos1 = junc_pos1,
      junc_pos2 = junc_pos2,
      junc_df = junc_df,
      transcript_db = toy_transcripts_gr
    )
  expect_true(nrow(transcripts_covering_junction) == 1)
  expect_true(ncol(transcripts_covering_junction) == 8)
  expect_true(transcripts_covering_junction$enst[1] == "ENST00000434301.1_1")
})
