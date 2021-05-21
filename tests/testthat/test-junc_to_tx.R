test_that("junc_to_tx works with custom toy data", {

  # example data
  transcripts <- GenomicRanges::GRangesList(list(
    tx1 = GenomicRanges::GRanges(
      c("1", "1", "1", "1"),
      IRanges::IRanges(
        c(5, 12, 18, 25),
        c(8, 14, 21, 26),
      ),
      strand = c("+", "+", "+", "+")
    ),
    tx2 = GenomicRanges::GRanges(
      c("1", "1"),
      IRanges::IRanges(
        c(5, 18),
        c(8, 21),
      ), strand = c("-", "-")
    ),
    tx3 = GenomicRanges::GRanges(
      c("1", "1"),
      IRanges::IRanges(
        c(5, 12),
        c(8, 14),
      ), strand = c("+", "+")
    )
  ))

  # junction examples
  junc_df = tibble::tribble(
    ~"chr", ~"pos1", ~"pos2", ~"strand", ~"comment",

    "1",  2, 12, "+", "intron-exon",
    "1",  7, 12, "+", "alt 5prime in exon, pos",
    "1",  8, 15, "-", "alt 5prime in intron, neg",
    "1",  8, 18, "+", "exon skipping, pos",
    "1",  8, 9, "+", "intron retention at 5prime, pos",
    "1",  11, 12, "+", "intron retention at 3prime, pos",
    "1",  8, 10, "+", "alt 3prim in intron, pos",
    "1",  8, 13, "-", "alt 5prim in exon, neg",
    "1",  6, 7, "+", "adjacent within exon, no effect!"
  ) %>%
    mutate(
      junc_id = str_c(chr, pos1, pos2, strand, sep = "_")
    )

  # junc_gr <- GenomicRanges::GRanges(
  #   junc_df$chr,
  #   IRanges::IRanges(junc_df$pos1, junc_df$pos2),
  #   strand = junc_df$strand
  #   )

  # Test with toy data ---------------------------------------------------------

  # single junction
  junc_idx <- 3
  chr = junc_df[junc_idx, ]$chr
  pos1 = junc_df[junc_idx, ]$pos1
  pos2 = junc_df[junc_idx, ]$pos2
  # strand = junc_df$strand[[2]]
  # transcript_ragnes <- range(transcripts)

  tx_df <- junc_to_tx(chr, pos1, pos2, transcripts)

  expect_true(nrow(tx_df) > 0)
  expect_true(ncol(tx_df) > 0)
  expect_true(length(tx_df$tx[[1]]) <= length(transcripts))

})

test_that("junc_to_tx works with provided toy data", {

  junc_id <- toy_junc_id[[1]]

  chr = str_split_fixed(junc_id, "_", 4)[1]
  pos1 = str_split_fixed(junc_id, "_", 4)[2] %>% as.integer()
  pos2 = str_split_fixed(junc_id, "_", 4)[3] %>% as.integer()

  tx_df <- junc_to_tx(chr, pos1, pos2, toy_transcripts)

  expect_true(nrow(tx_df) > 0)
  expect_true(ncol(tx_df) > 0)
  expect_true(length(tx_df$tx[[1]]) <= length(toy_transcripts))


})
