

test_that("spladder.transform works", {
  path <-  system.file("extdata", "", package = "splice2neo")
  dat.combined <- spladder_transform(path)
  expect_true(ncol(dat.combined) == 7)
  expect_true(nrow(dat.combined) == 26)
  expect_true(all(!is.na(dat.combined$junc_id)))

})

test_that("spladder.transform works correctly", {
  path <-  system.file("extdata", "", package = "splice2neo")
  dat.combined <- spladder_transform(path)

  expected_junc_ids <-
    c(
      "chr1:896180-896673:+",
      "chr1:896180-897009:+",
      "chr1:889462-891303:-",
      "chr1:889903-891303:-",
      "chr1:763229-764383:+",
      "chr1:763233-764383:+",
      "chr1:934993-935072:-",
      "chr1:934993-935246:-",
      "chr1:764484-776580:+",
      "chr1:776753-783034:+",
      "chr1:764484-783034:+",
      "chr1:1018367-1019295:-",
      "chr1:1019466-1019733:-",
      "chr1:1018367-1019733:-",
      "chr1:565845-565846:+",
      "chr1:565869-565870:+",
      "chr1:17368-17369:-",
      "chr1:17605-17606:-",
      "chr1:2323397-2324671:+",
      "chr1:2324794-2327223:+",
      "chr1:2323397-2324002:+",
      "chr1:2324670-2327223:+",
      "chr1:59133593-59134304:-",
      "chr1:59134354-59134656:-",
      "chr1:59133593-59134102:-",
      "chr1:59134238-59134656:-"
    )

  expect_true(all(dat.combined$junc_id == expected_junc_ids))

})


test_that("spladder_import works", {
  path <-  system.file("extdata", "", package = "splice2neo")
  dat.combined <- import_spladder(path)
  expect_true(length(dat.combined) == 5)
  expect_true("mutex_exons" %in% names(dat.combined))


})

test_that("spladder_import fails returns error message when wrong folder", {

  expect_error(import_spladder("."))

})

