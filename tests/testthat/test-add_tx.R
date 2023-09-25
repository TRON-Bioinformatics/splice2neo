test_that("add_tx works on toy example data", {


  junc_df <- tibble(
    junc_id = na.omit(toy_junc_id)
  )

  junc_tx <- add_tx(junc_df, toy_transcripts)

  expect_true(all(na.omit(toy_junc_id) %in% junc_tx$junc_id))
  expect_true(nrow(junc_tx) >= nrow(junc_df))
  expect_equal(typeof(junc_tx$tx_lst), "list")

})

test_that("add_tx works on manaul toy example data", {


  junc_df <- tibble(
    junc_id = c("chr2:152389996-152392205:-", "chr2:152389996-152390729:-",
                "chr2:152389955-152389956:-")
  )

  junc_tx <- add_tx(junc_df, toy_transcripts)

  expect_true(all(c("chr2:152389996-152392205:-", "chr2:152389996-152390729:-",
                    "chr2:152389955-152389956:-") %in% junc_tx$junc_id))
  expect_true(nrow(junc_tx) >= nrow(junc_df))
  expect_equal(typeof(junc_tx$tx_lst), "list")

})
