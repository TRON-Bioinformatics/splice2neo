

test_that("generate_combined_dataset works", {
  path <-  system.file("extdata", "", package = "splice2neo")
  spladder_juncs   <- spladder_transform(path)

  path <- system.file("extdata", "test_regtools_Aligned.out.sorted.bam.junc", package = "splice2neo")
  regtools_juncs <- regtools_transform(path)

  dat.combined <- generate_combined_dataset(list("spladder" = spladder_juncs,
                                                 "regtools" = regtools_juncs))
  expect_true(ncol(dat.combined) == 12)
  expect_true(nrow(dat.combined) == 31)
  expect_true(all(!is.na(dat.combined$junc_id)))

})


test_that("add_identified_in_RNA works", {

  # this is still not optimal as test data for mutations and rna junctions do not fit together
  # we can only test that function does not break
  path_to_spladder <- system.file("extdata", "", package = "splice2neo")
  path_to_regtools <- system.file("extdata", "test_regtools_Aligned.out.sorted.bam.junc", package = "splice2neo")
  dat.combined <- add_identified_in_RNA(mutation_juncs = toy_junc_df,
                                        path_to_spladder = path_to_spladder,
                                        path_to_regtools = path_to_regtools)


  expect_true(ncol(dat.combined) == 4)
  expect_true(nrow(dat.combined) == nrow(toy_junc_df))

})

