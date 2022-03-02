

test_that("leafcutter count transformation works", {
  count_file <- system.file("extdata", "test_perind.counts.gz", package = "splice2neo")
  leafcutter_counts <- import_leafcutter_counts(count_file)

  expect_true(nrow(leafcutter_counts) == 3)
  expect_true(ncol(leafcutter_counts) == 2)

  cts <- transform_leafcutter_counts(leafcutter_counts)

  expect_true(nrow(cts) == nrow(leafcutter_counts))
  expect_true(ncol(cts) == 8)


})




test_that("complete leafcutter import function works", {

  path <-  system.file("extdata", "", package = "splice2neo")
  df <- leafcutter_transform(path)

  expect_true(ncol(df) == 8)
  expect_true(nrow(df) == 3)
  expect_true(all(!is.na(df$junc_id)))

})

test_that("complete leafcutter import function works correctly", {

  path <-  system.file("extdata", "", package = "splice2neo")
  df <- leafcutter_transform(path)

  expected_junc_ids <- c("chr1:896180-897009:+",
                          "chr1:889462-889844:-",
                          "chr1:889462-891303:-")

  expect_true(all(df$junc_id == expected_junc_ids))

})
