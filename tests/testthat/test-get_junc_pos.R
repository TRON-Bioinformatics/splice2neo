
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


test_that("get_junc_pos works with negativ strand", {

  jx <- GenomicRanges::GRanges("1:8-15:-")
  tx <- GenomicRanges::GRangesList(
    GenomicRanges::GRanges(c("1:5-8:-", "1:18-21:-"))
  )
  #     0        1         2
  #     123456789012345678901234567890
  #alt      ====      =======
  #pos +    1234
  #pos -              7654321

  tx_alt <- modify_tx(tx, jx)
  pos <- get_junc_pos(tx_alt, jx)
  expect_equal(pos, 7)

  # convert strand to +
  BiocGenerics::strand(tx_alt) <- "+"
  BiocGenerics::strand(jx) <- "+"
  pos <- get_junc_pos(tx_alt, jx)
  expect_equal(pos, 4)

})


test_that("get_junc_pos works on toy example when junc is not on transcipt", {
  jx <- GenomicRanges::GRanges(rep("chr1:15", 3))
  tx <- GenomicRanges::GRangesList(
    tx1 = GenomicRanges::GRanges(c(
      "chr1:10-20",
      "chr1:50-60",
      "chr1:80-100")),
    tx2 = GenomicRanges::GRanges(c(
      "chr1:10-20",
      "chr1:50-60",
      "chr1:80-120")),
    tx3 = GenomicRanges::GRanges(c(
      "chr1:50-60",
      "chr1:80-100"))
  )

  pos <- get_junc_pos(tx, jx)

  expect_true(!is.na(pos[1]))
  expect_true(!is.na(pos[2]))
  expect_true(is.na(pos[3]))

})


test_that("get_junc_pos works when an empty range is in tx", {

  tx <- GenomicRanges::GRangesList(
      list(
        GenomicRanges::GRanges(), # empty range
        GenomicRanges::GRanges(c(
          "1:2-3:+",
          "1:5-8:+",
          "1:10-12:+"))
    )
  )

  jx <- GenomicRanges::GRanges(c(
    "1:3-10:+",  # any
    "1:3-11:+"  # any
  ))

  pos <- get_junc_pos(tx, jx)

  expect_true(is.na(pos[1]))
  expect_equal(length(tx), length(pos))

  # Multiple empty ranges ------------------------------------------------------
  tx <- GenomicRanges::GRangesList(
    list(
      GenomicRanges::GRanges(), # empty range
      GenomicRanges::GRanges(), # empty range
      GenomicRanges::GRanges(c(
        "1:2-3:+",
        "1:5-8:+",
        "1:10-12:+"))
    )
  )

  jx <- GenomicRanges::GRanges(c(
    "1:3-9:+",  # any
    "1:3-10:+",  # any
    "1:3-11:+"  # any
  ))

  pos <- get_junc_pos(tx, jx)

  expect_true(is.na(pos[1]))
  expect_equal(length(tx), length(pos))
})


test_that("get_intronretention_alt_pos works on positive strand", {

  tx <- GenomicRanges::GRangesList(
      GenomicRanges::GRanges(c(
        "1:2-10:+",
        "1:15-20:+",
        "1:25-30:+"))
    ,
    GenomicRanges::GRanges(c(
      "1:2-10:+",
      "1:15-20:+",
      "1:25-30:+")
      ),
    GenomicRanges::GRanges(c(
      "1:2-10:+",
      "1:15-20:+",
      "1:25-30:+")
    ))



  jx <- GenomicRanges::GRanges(c(
    "1:10-14:+",  # any
    "1:10-11:+" ,
    "1:14-15:+" # any
  ))

  tx_mod <- modify_tx(tx, jx)
  intron_retention <- ifelse(jx@ranges@width == 2, TRUE, FALSE)
  pos <- get_junc_pos(tx_mod, jx)
  other_pos <- get_intronretention_alt_pos(tx, tx_mod, jx, intron_retention)

  expect_true(other_pos[2] == 14)
  expect_true(other_pos[3] == 9)
  expect_true(other_pos[3] == pos[2])
  expect_true(other_pos[2]-1 == pos[3])

})

test_that("get_intronretention_alt_pos works on negative strand", {

  tx <- GenomicRanges::GRangesList(
    GenomicRanges::GRanges(c(
      "1:2-10:-",
      "1:15-20:-",
      "1:25-30:-"))
    ,
    GenomicRanges::GRanges(c(
      "1:2-10:-",
      "1:15-20:-",
      "1:25-30:-")
    ),
    GenomicRanges::GRanges(c(
      "1:2-10:-",
      "1:15-20:-",
      "1:25-30:-")
    ))



  jx <- GenomicRanges::GRanges(c(
    "1:10-14:-",  # any
    "1:10-11:-" ,
    "1:14-15:-" # any
  ))

  tx_mod <- modify_tx(tx, jx)
  intron_retention <- ifelse(jx@ranges@width == 2, TRUE, FALSE)
  pos <- get_junc_pos(tx_mod, jx)
  other_pos <- get_intronretention_alt_pos(tx, tx_mod, jx, intron_retention)

  expect_true(other_pos[2] == 12)
  expect_true(other_pos[3] == 17)
  expect_true(other_pos[2] == pos[3])
  expect_true(other_pos[3]-1 == pos[2])

})



