

test_that("leafcutter count transformation works", {
  count_file <- system.file("extdata", "test_perind.counts.gz", package = "splice2neo")
  leafcutter_counts <- import_leafcutter_counts(count_file)

  expect_true(nrow(leafcutter_counts) == 3)
  expect_true(ncol(leafcutter_counts) == 2)

  cts <- transform_leafcutter_counts(leafcutter_counts)

  expect_true(nrow(cts) == nrow(leafcutter_counts))
  expect_true(ncol(cts) == 6)


})


test_that("leafcutter bam.junc transformation works", {
  junc_file <- system.file("extdata", "test_Aligned.out.sorted.bam.junc", package = "splice2neo")
  leafcutter_bam_junc <- import_leafcutter_bam(junc_file)

  expect_true(ncol(leafcutter_bam_junc) == 12)
  expect_true(all(unique(leafcutter_bam_junc$strand) %in% c("+", "-")))

  jcs <- transform_leafcutter_bam(leafcutter_bam_junc)

  expect_true(ncol(jcs) == 2)
  expect_true(all(unique(jcs$strand) %in% c("+", "-")))

})




test_that("complete leafcutter import function works", {

  path <-  system.file("extdata", "", package = "splice2neo")
  df <- leafcutter_transform(path)

  expect_true(ncol(df) == 8)
  expect_true(nrow(df) == 3)
  expect_true(all(!is.na(df$junc_id)))

})
