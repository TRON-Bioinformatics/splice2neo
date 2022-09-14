test_that("parse.irfinder.txt works", {
  path <-  system.file("extdata", "", package = "splice2neo")
  dat.combined <- parse_irfinder_txt(path, warnings = FALSE, irratio = 0.01, cnn=FALSE)
  expect_true(ncol(dat.combined) == 12)
  expect_true(nrow(dat.combined) == 8)
  expect_true(all(c("overlapping_feature", "IRratio", "Warnings", "IntronDepth") %in% colnames(dat.combined)))

})

test_that("parse.irfinder.txt works correctly", {

  path <-  system.file("extdata", "", package = "splice2neo")
  dat.combined <- parse_irfinder_txt(path, warnings = FALSE, irratio = 0.01, cnn=FALSE)

  expected_junc_ids <- c("chr6:96538810-96538811:+",
                         "chr6:96540534-96540535:+",
                         "chr7:100673923-100673924:+",
                         "chr7:100676176-100676177:+",
                         "chr2:231791808-231791809:+",
                         "chr2:231794262-231794263:+",
                         "chr7:149547776-149547777:-",
                         "chr7:149549192-149549193:-")

  expect_true(all(dat.combined$junc_id == expected_junc_ids))
  expect_true(all(dat.combined$Warnings == '-'))
  expect_true(all(dat.combined$IRratio > 0.01))

})

test_that("filter.irfinder.txt works with warnings", {

  path <-  system.file("extdata", "", package = "splice2neo")
  dat_with_warnings <- parse_irfinder_txt(path, warnings = TRUE, irratio = 0.01, cnn=FALSE)
  expect_true(ncol(dat_with_warnings) == 12)
  expect_true(nrow(dat_with_warnings) == 12)
  expect_false(all(dat_with_warnings$Warnings == "-"))

})

test_that("parse.irfinder.txt fails and returns error message when IRFinder-IR-nondir.txt is missing", {
  path <-  "."
  expect_error(parse_irfinder_txt(path))

})

