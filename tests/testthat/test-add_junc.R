test_that("add_junc works with toy example data", {

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

  # Test with toy data ---------------------------------------------------------
  tx <- transcripts[c(1, 3)]
  jx <- GenomicRanges::GRanges(junc_df$chr,
                               IRanges::IRanges(junc_df$pos1, junc_df$pos2),
                                strand = junc_df$strand)[c(2, 2)]

  tx_alt <- add_junc(tx, jx)

  expect_equal(length(tx_alt), length(tx))
  expect_equal(end(tx_alt[[1]])[1], 7)

})

test_that("grl_update_end_at and grl_update_end_at works with toy example data", {

  tx <- GenomicRanges::GRangesList(list(
    GenomicRanges::GRanges(c("1:5-8:+", "1:12-14:+", "1:18-21:+", "1:25-26:+")),
    GenomicRanges::GRanges(c("1:5-8:+", "1:12-14:+"))
  ))


  at = c(1, 1)
  pos = c(7, 7)

  # modify end -----------------------------------------------------------------
  tx_alt_end <- grl_update_end_at(tx, at, pos)

  end_alt = tx_alt_end %>%
    as.list() %>%
    purrr::map2(at, magrittr::extract) %>%
    purrr::map_int(BiocGenerics::end)
  names(end_alt) <- NULL

  expect_equal(end_alt, pos)

  # modify start -----------------------------------------------------------------
  tx_alt_start <- grl_update_start_at(tx, at, pos)

  start_alt = tx_alt_start %>%
    as.list() %>%
    purrr::map2(at, magrittr::extract) %>%
    purrr::map_int(BiocGenerics::start)
  names(start_alt) <- NULL

  expect_equal(start_alt, pos)

})



