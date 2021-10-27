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
    dplyr::mutate(
      junc_id = stringr::str_c(chr, pos1, pos2, strand, sep = "_")
    )

  # Test with toy data ---------------------------------------------------------
  tx <- transcripts[c(1, 3)]
  jx <- GenomicRanges::GRanges(junc_df$chr,
                               IRanges::IRanges(junc_df$pos1, junc_df$pos2),
                                strand = junc_df$strand)[c(2, 2)]

  tx_alt <- add_junc(tx, jx)

  expect_equal(length(tx_alt), length(tx))
  expect_equal(S4Vectors::end(tx_alt[[1]])[1], 7)

})

test_that("add_junc work for exon exclusion junctions on toy data", {

  tx <- GenomicRanges::GRangesList(list(
    GenomicRanges::GRanges(c("1:2-3:+",
                             "1:5-6:+",
                             "1:10-15:+"))
    ))

  jx <- GenomicRanges::GRanges(c("1:3-10:+"))

  tx_alt <- add_junc(tx, jx)

  expect_equal(length(tx_alt), length(tx))
  expect_equal(tx_alt[[1]], tx[[1]][c(1, 3)])

})


test_that("add_junc work for exitorn junction on toy data", {

  tx <- GenomicRanges::GRangesList(list(
    GenomicRanges::GRanges(c("1:2-3:+",
                             "1:10-30:+",
                             "1:40-50:+"))
  ) %>% rep(3))

  jx <- GenomicRanges::GRanges(c(
    "1:3-10:+",  # canonical
    "1:15-20:+", # exitron in second exon
    "1:30-40:+"  # canonical
  ))

  tx_alt <- add_junc(tx, jx)

  expect_equal(length(tx_alt[[2]]), length(tx[[2]]) + 1)
  expect_equal(tx_alt[[1]], tx[[1]])
  expect_equal(tx_alt[[3]], tx[[3]])
  expect_equal(length(tx_alt), length(tx))

})


test_that("add_junc_pos works on toy data", {

  tx <- GenomicRanges::GRangesList(
    c(
      list(
        GenomicRanges::GRanges(c("1:2-3:+",
                             "1:10-30:+",
                             "1:40-50:+"))
      ) %>% rep(3),
      GenomicRanges::GRanges(c("1:2-3:-",
                               "1:5-8:-",
                               "1:10-12:-"))
    )
  )

  jx <- GenomicRanges::GRanges(c(
    "1:3-10:+",  # canonical
    "1:15-20:+", # exitron in second exon
    "1:30-40:+", # canonical
    "1:7-10:-"  # canonical neg. strand
  ))

  pos <- add_junc_pos(tx, jx)

  expect_equal(length(tx), length(pos))

})




