test_that("canonical_junctions works on test gtf", {


  gtf_file <- system.file("extdata","GTF_files","Aedes_aegypti.partial.gtf",
                          package="GenomicFeatures")

  tx <- parse_gtf(gtf_file, format = "gtf")

  jx <- canonical_junctions(tx)

  expect_equal(typeof(jx), "character")
  expect_true(length(jx) >= 1)
  expect_equal(length(jx), length(unique(jx)))

})
