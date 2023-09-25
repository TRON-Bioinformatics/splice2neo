

test_that("generate_combined_dataset works", {
  path <-  system.file("extdata", "", package = "splice2neo")
  spladder_juncs <- spladder_transform(path)
  leafcutter_juncs <- leafcutter_transform(path)

  dat.combined <- generate_combined_dataset(spladder_juncs,
                                            leafcutter_juncs)
  expect_true(ncol(dat.combined) == 9)
  expect_true(nrow(dat.combined) == 27)
  expect_true(all(!is.na(dat.combined$junc_id)))
})


test_that("add_identified_in_RNA works", {

  # this is still not optimal as test data for mutations and rna junctions do not fit together
  # we can only test that function does not break
  path <-  system.file("extdata", "", package = "splice2neo")
  dat.combined <- add_identified_in_RNA(mutation_juncs = toy_junc_df, path_to_spladder = path, path_to_leafcutter = path)


  expect_true(ncol(dat.combined) == 4)
  expect_true(nrow(dat.combined) == nrow(toy_junc_df))

})

