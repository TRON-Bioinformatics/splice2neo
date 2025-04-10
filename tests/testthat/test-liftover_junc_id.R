test_that("liftover_junc_id works on toy example data", {

  chain_file = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
  stopifnot(chain_file != "")

  junc_df <- toy_junc_df

  junc_df_lifted <- liftover_junc_id(junc_df, chain_file)

  expect_equal(nrow(junc_df_lifted), nrow(junc_df))
  expect_true(all(c("liftover_successful", "liftover_unique", "junc_id_lifted") %in% names(junc_df_lifted)))


})



test_that("liftover_junc_id works with non-unique mappings", {

  chain_file = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
  stopifnot(chain_file != "")

  junc_df <- toy_junc_df %>%
    head(3) %>%
    add_row(junc_id = "chr2:10000000-10050000:-") %>%
    add_row(junc_id = "chr7:123-456:+")

  junc_df_lifted <- liftover_junc_id(junc_df, chain_file)

  expect_equal(nrow(junc_df_lifted), nrow(junc_df))
  expect_false(all(junc_df_lifted$liftover_unique))
  expect_false(all(junc_df_lifted$liftover_successful))
  expect_true(all(junc_df_lifted$liftover_unique[1:3]))
  expect_true(any(is.na(junc_df_lifted$junc_id_lifted)))

})


test_that("liftover_junc_id for deleted junction position", {

  chain_file = system.file(package="splice2neo", "extdata", "test.chain")
  stopifnot(chain_file != "")

  junc_df <- dplyr::tibble(
    junc_id = "chr6:35233301-35233342"
  )

  junc_df_lifted <- liftover_junc_id(junc_df, chain_file)

  expect_equal(nrow(junc_df_lifted), nrow(junc_df))
  expect_false(all(junc_df_lifted$liftover_successful))
  expect_true(all(junc_df_lifted$liftover_unique))

})


test_that("liftover_junc_id for deletion located at the beginning of the junction", {

  chain_file = system.file(package="splice2neo", "extdata", "test.chain")
  stopifnot(chain_file != "")

  junc_df <- dplyr::tibble(
    junc_id = "chr6:30704795-30705795"
  )

  junc_df_lifted <- liftover_junc_id(junc_df, chain_file)

  expect_equal(nrow(junc_df_lifted), nrow(junc_df))
  expect_true(all(junc_df_lifted$liftover_successful))
  expect_false(all(junc_df_lifted$liftover_unique))

})

