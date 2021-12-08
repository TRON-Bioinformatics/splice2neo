test_that("bed_to_junc works on rtracklayer test file", {

  test_bed <- system.file("tests","test.bed", package = "rtracklayer")
  junc_ids <- bed_to_junc(test_bed)

  expect_true(length(junc_ids) == 5)

})


test_that("bed_to_junc works on rtracklayer test file with intron", {

  test_bed <- system.file("tests","test.bed", package = "rtracklayer")
  junc_ids <- bed_to_junc(test_bed, type = "intron")

  expect_true(length(junc_ids) == 5)

})
