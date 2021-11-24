test_that("transform_for_requant works", {
  cts_id <- c("0cb012d5e", "0cb012d5e", "1276125bbd07")
  cts_seq <- c("GCAAATATGGATTAAGCCGC", "GCAAATATGGATTAAGCCGC", "CCGCCAAAGCATTATGGAG")
  cts_junc_pos <- c(7, 7, 7)
  df <- tibble(cts_id , cts_seq, cts_junc_pos)

  df_easy <- transform_for_requant(df)
  expect_equal(nrow(df_easy), 2)
  expect_equal(ncol(df_easy), 3)


  cts_id <- c("0cb012d5e", "0cb012d5e", "1276125bbd07")
  cts_seq <- c("GCAAATATGGATTAAGCCGC", NA, "CCGCCAAAGCATTATGGAG")
  cts_junc_pos <- c(7, NA, 7)
  df <- tibble(cts_id , cts_seq, cts_junc_pos)

  df_easy <- transform_for_requant(df)
  expect_equal(nrow(df_easy), 2)
  expect_equal(ncol(df_easy), 3)
})

