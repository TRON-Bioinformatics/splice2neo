

test_that("regtools.transform works", {
  path <- system.file("extdata", "test_regtools_Aligned.out.sorted.bam.junc", package = "splice2neo")
  df <- regtools_transform(path)
  expect_true(ncol(df) == 9)
  expect_true(nrow(df) == 5)
  expect_true(all(!is.na(df$junc_id)))

})

test_that("regtools.transform works correctly", {
  path <-  system.file("extdata", "test_regtools_Aligned.out.sorted.bam.junc", package = "splice2neo")
  df <- regtools_transform(path)
  expected_junc_ids <-
    c(
      "chr1:14829-14970:-",
      "chr1:14843-185491:-",
      "chr1:17055-17606:-",
      "chr1:17368-17606:-",
      "chr1:17368-17526:-"
    )
  expected_classification <-
    c(
      "novel_acceptor_novel_donor",
      "A3SS",
      "cassette_exon",
      "canonical_junction",
      "A5SS"
    )
  expect_true(all(df$junc_id == expected_junc_ids))
  expect_true(all(df$class == expected_classification))
})


test_that("import_regtools_junc fails and returns error message when file doesn't exists", {

  expect_error(import_regtools_junc("./test_regtools_Aligned.out.sorted.bam.junc"))

})
