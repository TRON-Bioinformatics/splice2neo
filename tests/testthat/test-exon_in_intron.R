test_that("exon_in_intron works in case there is a exon of another transcripts in intron", {

  toy_tx <- GenomicRanges::GRangesList(list(
    tx1 = GenomicRanges::GRanges(
      c("1", "1", "1"),
      IRanges::IRanges(
        c(2, 15, 27),
        c(8, 23, 35),
      ),
      strand = c("+", "+", "+")
    ),
    tx1a = GenomicRanges::GRanges(
      c("1" ),
      IRanges::IRanges(
        c(3),
        c(10),
      ),
      strand = c("+")
    ),
    tx2 = GenomicRanges::GRanges(
      c("1", "1", "1"),
      IRanges::IRanges(
        c(2, 15, 27),
        c(8, 23, 35),
      ),
      strand = c("-", "-", "-")
    ),
    tx2a = GenomicRanges::GRanges(
      c("1"),
      IRanges::IRanges(
        c(10),
        c(13),
      ),
      strand = c("-")
    )
  ))

  df <- tibble(
    junc_id = c("1:8-9:+", "1:14-15:+",  "1:8-9:-", "1:14-15:-"),
    tx_id = c("tx1", "tx1", "tx2", "tx2")
  )

  df1 <- df %>%
    exon_in_intron( transcripts = toy_tx)

  expect_equal(unique(df1$exon_free), FALSE)
  expect_equal(nrow(df), nrow(df1))
})


test_that("exon_in_intron works in case there is NO exon of another transcripts in intron", {

  toy_tx <- GenomicRanges::GRangesList(list(
    tx1 = GenomicRanges::GRanges(
      c("1", "1", "1"),
      IRanges::IRanges(
        c(2, 15, 27),
        c(8, 23, 35),
      ),
      strand = c("+", "+", "+")
    ),
    tx1a = GenomicRanges::GRanges(
      c("1" ),
      IRanges::IRanges(
        c(3),
        c(10),
      ),
      strand = c("+")
    ),
    tx2 = GenomicRanges::GRanges(
      c("1", "1", "1"),
      IRanges::IRanges(
        c(2, 15, 27),
        c(8, 23, 35),
      ),
      strand = c("-", "-", "-")
    ),
    tx2a = GenomicRanges::GRanges(
      c("1"),
      IRanges::IRanges(
        c(10),
        c(13),
      ),
      strand = c("-")
    )
  ))

  df <- tibble(
    junc_id = c("1:23-24:+", "1:26-27:+",  "1:23-24:-", "1:26-27:-"),
    tx_id = c("tx1", "tx1", "tx2", "tx2")
  )

  df1 <- df %>%
    exon_in_intron(transcripts = toy_tx)

  expect_equal(unique(df1$exon_free), TRUE)
  expect_equal(nrow(df), nrow(df1))
})


test_that("exon_in_intron works for non IRs", {

  toy_tx <- GenomicRanges::GRangesList(list(
    tx1 = GenomicRanges::GRanges(
      c("1", "1", "1"),
      IRanges::IRanges(
        c(2, 15, 27),
        c(8, 23, 35),
      ),
      strand = c("+", "+", "+")
    ),
    tx1a = GenomicRanges::GRanges(
      c("1" ),
      IRanges::IRanges(
        c(3),
        c(10),
      ),
      strand = c("+")
    ),
    tx2 = GenomicRanges::GRanges(
      c("1", "1", "1"),
      IRanges::IRanges(
        c(2, 15, 27),
        c(8, 23, 35),
      ),
      strand = c("-", "-", "-")
    ),
    tx2a = GenomicRanges::GRanges(
      c("1"),
      IRanges::IRanges(
        c(10),
        c(13),
      ),
      strand = c("-")
    )
  ))

  df <- tibble(
    junc_id = c("1:8-9:+", "1:23-24:+", "1:8-13:+", "1:8-13:-"),
    tx_id = c("tx1", "tx1", "tx1", "tx2")
  )

  df1 <- df %>%
    exon_in_intron( transcripts = toy_tx)

  expect_true(is.na(unique(df1$exon_free[3:4])))
  expect_equal(nrow(df), nrow(df1))
})

test_that("exon_in_intron does not fail for junctions at the end of a transcript ", {

  toy_tx <- GenomicRanges::GRangesList(list(
    tx1 = GenomicRanges::GRanges(
      c("1", "1", "1"),
      IRanges::IRanges(
        c(2, 15, 27),
        c(8, 23, 35),
      ),
      strand = c("+", "+", "+")
    ),
    tx1a = GenomicRanges::GRanges(
      c("1" ),
      IRanges::IRanges(
        c(3),
        c(10),
      ),
      strand = c("+")
    ),
    tx2 = GenomicRanges::GRanges(
      c("1", "1", "1"),
      IRanges::IRanges(
        c(2, 15, 27),
        c(8, 23, 35),
      ),
      strand = c("-", "-", "-")
    ),
    tx2a = GenomicRanges::GRanges(
      c("1"),
      IRanges::IRanges(
        c(10),
        c(13),
      ),
      strand = c("-")
    )
  ))

  df <- tibble(
    junc_id = c("1:35-36:+", "1:35-36:-"),
    tx_id = c("tx1",  "tx2")
  )

  df1 <- df %>%
    exon_in_intron(transcripts = toy_tx)

  expect_true(is.na(unique(df1$exon_free)))
  expect_equal(nrow(df), nrow(df1))
})

