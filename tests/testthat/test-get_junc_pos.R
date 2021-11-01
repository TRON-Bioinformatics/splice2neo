
test_that("get_junc_pos works on toy data", {

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

  pos <- get_junc_pos(tx, jx)

  expect_equal(length(tx), length(pos))

})


