test_that("format_spliceai works on SpliceAI example file", {

  spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")
  df_raw <- parse_spliceai(spliceai_file)
  df <- format_spliceai(df_raw)

  expect_true(nrow(df) >= 10)
  expect_true(!"gene_id" %in% names(df))

})

test_that("format_spliceai works on SpliceAI example file with gene_table", {

  gene_table <- tibble::tibble(
    gene_name = c("NEB", "TTN", "RYR1", "COL6A1" ,"AIFM1" ),
    gene_id = c("id1", "id2", "id3", "id4", "id5")
  )

  spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")
  df_raw <- parse_spliceai(spliceai_file)
  df <- format_spliceai(df_raw, gene_table = gene_table)

  expect_true(nrow(df) >= 10)
  expect_true("gene_id" %in% names(df))
  expect_true(all(!is.na(df$gene_id)))

})

test_that("format_spliceai_thresh works on SpliceAI example file", {

  spliceai_file <- system.file("extdata", "spliceai_thresh_output.vcf", package = "splice2neo")
  df_raw <- parse_spliceai_thresh(spliceai_file)
  df <- format_spliceai_thresh(df_raw)

  expect_true(nrow(df) >= 10)

})
