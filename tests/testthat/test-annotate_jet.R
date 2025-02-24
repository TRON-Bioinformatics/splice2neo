test_that("Annotate of potential JETs works", {

  rmsk <- readr::read_tsv(system.file("extdata", "rmsk_hg19_subset.tsv.gz", package = "splice2neo"))
  rmsk <- GenomicRanges::makeGRangesFromDataFrame(rmsl)

  junc_df <- tibble::tibble(
    junc_id = c("chr2:152407458-152408252:-")
  )

  df_jet <- annotate_potential_jet(junc_df, rmsk)
  
  expect_true(length(df_jet) ==  1)
  expect_true(df_jet$potential_jet)
  expect_true(df_jet$left_side_retroelement == "AluJb")
  expect_true(df_jet$left_side_retroelement_class == "SINE")
  expect_true(is.na(df_jet$right_side_retroelement_class))
  expect_true(is.na(df_jet$right_side_retroelement_class))


})