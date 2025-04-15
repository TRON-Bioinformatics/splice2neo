test_that("choose_tx works on toy example data and provides expected results ", {


  junc_df <- tibble(
    junc_id = na.omit(toy_junc_id)
  )

  df <- add_tx(junc_df, toy_transcripts)

  df_filtered <- choose_tx(df)

  # works?
  expect_true(all(na.omit(toy_junc_id) %in% df_filtered$junc_id))
  expect_true(nrow(df) >= nrow(df_filtered))
  expect_true(all(df_filtered$putative_event_type %in% c("ASS", "ES", "exitron", "IR", "ref junction")))


  # provides expected results?
  df_filtered1 <- df_filtered %>%
    tidyr::separate(junc_id, into = c("chr_junc","start-end","strand_junc" ), sep = ":", remove = FALSE) %>%
    tidyr::separate(`start-end`, into = c("start_junc", "end_junc" ), sep = "-") %>%
    mutate(
      start_junc = as.numeric(start_junc),
      end_junc = as.numeric(end_junc)
    ) %>%
    mutate(intron_retention = ifelse(abs(start_junc - end_junc) == 1, TRUE, FALSE))

  expect_true(all(df_filtered1 %>% filter(intron_retention) %>% pull(putative_event_type) == "IR"))
  expect_true(all(df_filtered1 %>% filter(junc_id == "chr2:152389996-152390729:-") %>% pull(putative_event_type) == "ASS"))


})

test_that("choose_tx works for exon skippings", {


  junc_df <- tibble(
    # ES, skipped exon is not occuring in all transcripts
    junc_id = c("chr2:152388410-152390729:-", "chr2:157332719-157352556:+" )
  )

  df <- add_tx(junc_df, toy_transcripts)

  df_filtered <- choose_tx(df)

  expect_true(all(
    df_filtered %>% filter(junc_id == "chr2:157332719-157352556:+") %>% pull(putative_event_type) == c(
      "ref junction",
      "ref junction",
      "ref junction",
      "ref junction",
      "ref junction",
      "ES",
      "ref junction",
      "ref junction",
      "ref junction",
      "ES"
    )
  ))

  expect_true(all(
    df_filtered %>% filter(junc_id == "chr2:152388410-152390729:-") %>% pull(putative_event_type) == c(
      "ref junction",
      "ref junction",
      "ref junction",
      "ES",
      "ref junction",
      "ref junction"
    )
  ))



})

test_that("choose_tx annotates correctly the distance to canonical splice sites on the plus strand ", {


  # Fake multi exon transcript on plus strand
  toy_tx <- GenomicRanges::GRanges(
    c("chr1:1-3:+", "chr1:8-10:+", "chr1:15-17:+", "chr1:22-26:+", "chr1:28-35:+")
  )
  toy_tx <- GenomicRanges::GRangesList(toy_tx)
  names(toy_tx) <- "E1"


  toy_jx <- tibble(junc_id = c("chr1:2-8:+",   # ASS with donor in exon
                               "chr1:5-8:+",   # ASS with donor in intron
                               "chr1:10-14:+", # ASS with acceptor in intron
                               "chr1:10-16:+", # ASS with acceptor in exon
                               "chr1:17-28:+", # ES
                               "chr1:30-33:+", # Exitron
                               "chr1:10-11:+", # IR like junction
                               "chr1:10-24:+", # Complex - ES with ASS at acceptor end
                               "chr1:3-8:+"),  #
                   tx_lst = toy_tx %>% as.list(), tx_id = "E1")

  toy_jx <- toy_jx %>% choose_tx()

  # Works with test data?
  expect_true(all(toy_jx$putative_event_type %in% c("ASS", "ASS", "ASS", "ASS", "ES", "exitron", "IR", "ASS", "ref junction")))
  expect_true(all(toy_jx$distance_to_next_canonical_donor == c(1, 2, 0, 0, 0, 2, 0, 0, 0)))
  expect_true(all(toy_jx$distance_to_next_canonical_acceptor == c(0, 0, 1, 1, 0, 2, 0, 2, 0)))

})

test_that("choose_tx annotates correctly the distance to canonical splice sites on the minus strand ", {


  # Fake multi exon transcript on minus strand
  toy_tx <- GenomicRanges::GRanges(
    c("chr1:1-3:-", "chr1:8-10:-", "chr1:15-17:-", "chr1:22-26:-", "chr1:28-35:-")
  )
  toy_tx <- GenomicRanges::GRangesList(toy_tx)
  names(toy_tx) <- "E2"


  toy_jx <- tibble(junc_id = c("chr1:2-8:-",   # ASS with acceptor in exon
                               "chr1:5-8:-",   # ASS with acceptor in intron
                               "chr1:10-14:-", # ASS with donor in intron
                               "chr1:10-16:-", # ASS with donor in exon
                               "chr1:17-28:-", # ES
                               "chr1:30-34:-", # Exitron
                               "chr1:10-11:-", # IR like junction
                               "chr1:10-24:-"),# Complex - ES with ASS at donor end
                   tx_lst = toy_tx %>% as.list(), tx_id = "E2")

  toy_jx <- toy_jx %>% choose_tx()

  # Works with test data?
  expect_true(all(toy_jx$putative_event_type %in% c("ASS", "ASS", "ASS", "ASS", "ES", "exitron", "IR", "ASS")))
  expect_true(all(toy_jx$distance_to_next_canonical_donor == c(0, 0, 1, 1, 0, 1, 0, 2)))
  expect_true(all(toy_jx$distance_to_next_canonical_acceptor == c(1, 2, 0, 0, 0, 2, 0, 0)))

})


