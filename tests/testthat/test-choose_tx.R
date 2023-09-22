test_that("choose_tx works on toy example data and provides expected results ", {


  junc_df <- tibble(
    junc_id = na.omit(toy_junc_id)
  )

  df <- add_tx(junc_df, toy_transcripts)

  df_filtered <- choose_tx(df)

  # works?
  expect_true(all(na.omit(toy_junc_id) %in% df_filtered$junc_id))
  expect_true(nrow(df) >= nrow(df_filtered))
  expect_true(all(df_filtered$putative_event_type %in% c("ASS", "ES", "exitron", "IR", "WT junction")))


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
      "WT junction",
      "WT junction",
      "WT junction",
      "WT junction",
      "WT junction",
      "ES",
      "WT junction",
      "WT junction",
      "WT junction",
      "ES"
    )
  ))

  expect_true(all(
    df_filtered %>% filter(junc_id == "chr2:152388410-152390729:-") %>% pull(putative_event_type) == c(
      "WT junction",
      "WT junction",
      "WT junction",
      "ES",
      "WT junction",
      "WT junction"
    )
  ))



})

