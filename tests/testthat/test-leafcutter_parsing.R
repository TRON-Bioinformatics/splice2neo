

test_that("leafcutter count transformation works", {

  cts <- transform_leafcutter_counts(leafcutter_counts)

  expect_true(nrow(cts) == nrow(leafcutter_counts))
  expect_true("chr1_21877890_21878165" == cts$fake_junc_id[1])


})


test_that("leafcutter bam.junc transformation works", {
  jcs <- transform_leafcutter_bam(leafcutter_bam_junc)

  expect_true(ncol(jcs) == 2)
  expect_true(jcs$strand[1] == "-")

})


test_that("general leafcutter transformation works", {
  jcs <- transform_leafcutter_bam(leafcutter_bam_junc)
  cts <- transform_leafcutter_counts(leafcutter_counts)
  dt <- left_join(cts, jcs, by = "fake_junc_id")
  dat_out <- leafcutter_transform_format(dt)

  expect_true(ncol(dat_out) == 8)
  expect_true(nrow(dat_out) == 2)
  expect_true(dat_out$junc_id[1] == "chr1:21877890-21878165:-")

})
