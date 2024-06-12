test_that("parse_gtf works with GFF file", {

  gff_file <- system.file("extdata","GFF3_files","a.gff3",package="splice2neo")

  grl <- parse_gtf(gff_file, format = "gff3")

  expect_s4_class(grl, "CompressedGRangesList")
})

test_that("parse_gtf works with GTF file", {

  gtf_file <- system.file("extdata","GTF_files","Aedes_aegypti.partial.gtf",
                         package="splice2neo")

  grl <- parse_gtf(gtf_file, format = "gtf")

  expect_s4_class(grl, "CompressedGRangesList")
})

