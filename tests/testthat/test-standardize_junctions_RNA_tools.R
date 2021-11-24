test_that("generate_combined_dataset works", {
  dat.combined <- generate_combined_dataset(spladder_juncs,
                                            leafcutter_juncs)
  expect_true(ncol(dat.combined) == 10)
})


test_that("spladder.transform.format works", {
  dat.combined <- spladder_transform_format(spladder_output)
  expect_true(ncol(dat.combined) == 8)
})