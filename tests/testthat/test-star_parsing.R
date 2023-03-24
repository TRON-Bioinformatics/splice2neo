test_that("parse.star.sj works", {
  path <-  system.file("extdata", "test_star_SJ.out.tab", package = "splice2neo")
  dat.combined <- parse_star_sj(path)
  expect_true(ncol(dat.combined) == 10)
  expect_true(nrow(dat.combined) == 7)
  expect_true(all(!is.na(dat.combined$junc_id)))

})

test_that("parse.star.sj works correctly", {
  path <-  system.file("extdata", "test_star_SJ.out.tab", package = "splice2neo")
  dat.combined <- parse_star_sj(path)
  expected_junc_ids <-
    c(
      "chr1:12227-12613:+",
      "chr1:12697-13221:+",
      "chr1:12721-13221:+",
      "chr1:15038-15796:-",
      "chr1:15947-16607:-",
      "chrX:156021563-156021694:+",
      "chrX:156021792-156021999:+"
    )
  exepected_unknown_strand <-
    c("chr1:15940-16616:*",
      "chrX:156022145-156022314:*"
    )
  expect_true(all(dat.combined$junc_id == expected_junc_ids))
  expect_true(all(! exepected_unknown_strand %in% expected_junc_ids))
})

test_that("parse.star.sj returns parses raw RNA-seq support", {
  unique_reads <- c(0, 0, 0, 0, 0, 0, 14)
  multi_reads <- c(13, 7, 3, 61, 8, 9, 132)
  path <-  system.file("extdata", "test_star_SJ.out.tab", package = "splice2neo")
  dat.combined <- parse_star_sj(path)
  expect_true(all(dat.combined$uniquely_mapping_reads == unique_reads))
  expect_true(all(dat.combined$multi_mapping_reads == multi_reads))

})

test_that("import_star_sj returns error message when file doesn't exists", {

  expect_error(import_star_sj("./test_star_SJ.out.tab"))

})
