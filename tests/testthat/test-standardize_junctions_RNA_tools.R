test_that("generate_combined_dataset works", {
  dat.combined <- generate_combined_dataset(spladder_juncs,
                                            leafcutter_juncs, canonical_juncs)
  expect_true(ncol(dat.combined) == 12)
})


test_that("spladder.transform.format works", {
  dat.combined <- spladder.transform.format(spladder_output)
  expect_true(ncol(dat.combined) == 8)
})
