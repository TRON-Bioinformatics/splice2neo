test_that("parse.irfinder.txt works", {
  path <-  system.file("extdata", "", package = "splice2neo")
  dat.combined <- parse_irfinder_txt(path, warnings = FALSE, irratio = 0.01, cnn=FALSE)
  expect_true(ncol(dat.combined) == 12)
  expect_true(nrow(dat.combined) == 8)
  expect_true(all(c("overlapping_feature", "IRratio", "Warnings", "IntronDepth") %in% colnames(dat.combined)))

})

test_that("parse.irfinder.txt works with CNN input", {
  path <-  system.file("extdata", "", package = "splice2neo")
  dat.combined <- parse_irfinder_txt(path, warnings = FALSE, irratio = 0.01, cnn=TRUE)
  expect_true(ncol(dat.combined) == 12)
  expect_true(nrow(dat.combined) == 18)
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

test_that("parse.irfinder.txt works correctly with CNN input", {

  path <-  system.file("extdata", "", package = "splice2neo")
  dat.combined <- parse_irfinder_txt(path, warnings = FALSE, irratio = 0.01, cnn=TRUE)

  expected_junc_ids <- c("chr1:961552-961553:+",
                         "chr1:961628-961629:+",
                         "chr1:961750-961751:+",
                         "chr1:961825-961826:+",
                         "chr1:962047-962048:+",
                         "chr1:962354-962355:+",
                         "chr1:962471-962472:+",
                         "chr1:962703-962704:+",
                         "chr1:963253-963254:+",
                         "chr1:963336-963337:+",
                         "chr1:963504-963505:+",
                         "chr1:963919-963920:+",
                         "chr1:964008-964009:+",
                         "chr1:964106-964107:+",
                         "chr1:970758-970759:+",
                         "chr1:970878-970879:+",
                         "chr1:999613-999614:-",
                         "chr1:999691-999692:-")
  expect_true(all(dat.combined$junc_id == expected_junc_ids))
  expect_true(all(dat.combined$Warnings == '-'))
  expect_true(all(dat.combined$IRratio > 0.05))

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
  expect_error(parse_irfinder_txt(path, cnn=TRUE))


})

