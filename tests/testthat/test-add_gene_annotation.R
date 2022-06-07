test_that("add_gene_name works", {

  junc_df_annotated <- add_gene_name(toy_junc_df)

  expect_equal(nrow(junc_df_annotated), nrow(toy_junc_df))
  expect_true(ncol(junc_df_annotated) == ncol(toy_junc_df) + 1)
  expect_true(all(!is.na(junc_df_annotated$gene_name)))

})


test_that("add_gene_name works for coordinates with multiple genes", {

  df_test <- tibble(junc_id = "chr1:145497392-145497393:+")

  junc_df_annotated <- add_gene_name(df_test)

  expect_equal(nrow(junc_df_annotated), nrow(df_test))
  expect_true(ncol(junc_df_annotated) == ncol(df_test) + 1)
  expect_true(str_count(junc_df_annotated$gene_name,",") > 1)

})
