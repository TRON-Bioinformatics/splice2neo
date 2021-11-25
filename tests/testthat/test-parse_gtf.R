test_that("parse_gtf works with GFF file", {

  gff_file <- system.file("extdata","GFF3_files","a.gff3",package="GenomicFeatures")

  grl <- parse_gtf(gff_file)

  expect_s4_class(grl, "CompressedGRangesList")
})

test_that("parse_gtf works with GTF file", {

  gtf_file <- system.file("extdata","GTF_files","Aedes_aegypti.partial.gtf",
                         package="GenomicFeatures")

  grl <- parse_gtf(gtf_file)

  expect_s4_class(grl, "CompressedGRangesList")
})

