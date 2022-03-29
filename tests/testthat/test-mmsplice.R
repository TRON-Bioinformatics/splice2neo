test_that("parse_mmsplice works with example file", {

  mmsplice_file <- system.file("extdata", "mmsplice_pred.csv", package = "splice2neo")
  df <- parse_mmsplice(mmsplice_file)

  expect_true(nrow(df) >= 10)

})


test_that("get_exon_skipping_junction() works with toy example data", {

  mmsplice_file <- system.file("extdata", "mmsplice_pred.csv", package = "splice2neo")
  mmsplice_df <- parse_mmsplice(mmsplice_file)

  mmsplice_df <- mmsplice_df %>%
    dplyr::filter(transcript_id %in% names(toy_transcripts))

  junc_annot <- annotate_mmsplice(mmsplice_df, toy_transcripts)

  expect_true(nrow(junc_annot) >= 10)


})

test_that("get_exon_skipping_junction() and get_exon_sincljusion_junction work with toy data", {


  #=============================================================================
  #          1         2         3         4
  # 1234567890123456789012345678901234567890
  # -#######------#########---#########-----  tx1 (positive strand)
  #               xxxxxxxxx                   exon skipping
  #        [==================]               exon skipping junction
  #               iiiiiiiii                   exon inclusion
  #        [======]                           exon inclusion junction 1
  #                       [===]               exon inclusion junction 2
  # 1234567890123456789012345678901234567890
  #          1         2         3         4
  # -#######------#########---#########-----  tx2 (negative strand)
  #               xxxxxxxxx                   exon skipping
  #        [==================]               exon skipping junction
  #               iiiiiiiii                   exon inclusion
  #        [======]                           exon inclusion junction 1
  #                       [===]               exon inclusion junction 2
  #=============================================================================


  toy_tx <- GenomicRanges::GRangesList(list(
    tx1 = GenomicRanges::GRanges(
      c("1", "1", "1"),
      IRanges::IRanges(
        c(2, 15, 27),
        c(8, 23, 35),
      ),
      strand = c("+", "+", "+"),
      exon_name = c("ex1", "ex2", "ex3")
    ),
    tx2 = GenomicRanges::GRanges(
      c("1", "1", "1"),
      IRanges::IRanges(
        rev(c(2, 15, 27)),
        rev(c(8, 23, 35)),
      ),
      strand = c("-", "-", "-"),
      exon_name = rev(c("ex1", "ex2", "ex3"))
    )
  ))

  expected_pos = c("1:8-27:+",
                   "1:8-15:+",
                   "1:23-27:+"
                    )
  expected_neg = c("1:8-27:-",
                   "1:8-15:-",
                   "1:23-27:-"
                    )

  # run function and combine observed results
  obs_pos <- c(
    get_exon_skipping_junction("ex2", "tx1", toy_tx),
    get_exon_inclusion_junction("ex2", "tx1", toy_tx)
  ) %>% unlist()

  # run function and combine observed results for negative transcript
  obs_neg <- c(
    get_exon_skipping_junction("ex2", "tx2", toy_tx),
    get_exon_inclusion_junction("ex2", "tx2", toy_tx)
  ) %>% unlist()


  expect_equal(obs_pos, expected_pos)
  expect_equal(obs_neg, expected_neg)

  # test first and last exons which should lead to NA  -------------------------
  l1 <- get_exon_skipping_junction("ex1", "tx1", toy_tx)
  expect_equal(l1[[1]], NA_character_)

  l2 <- get_exon_inclusion_junction("ex1", "tx1", toy_tx)
  expect_equal(l2[[1]][1], NA_character_)
  expect_true(!is.na(l2[[1]][2]))

  # test vectorized use --------------------------------------------------------
  l1 <- get_exon_skipping_junction(c("ex2", "ex2"), c("tx1", "tx2"), toy_tx)
  l2 <- get_exon_inclusion_junction(c("ex2", "ex2"), c("tx1", "tx2"), toy_tx)
  expect_equal(length(l1), 2)
  expect_equal(length(l2), 2)

})
