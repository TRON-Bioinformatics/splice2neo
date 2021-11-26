test_that("add_context_seq works on toy example data", {

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  junc_df <- tibble::tibble(
    junc_id = toy_junc_id[c(1, 6, 10)]
  ) %>%
    add_tx(toy_transcripts)

  cts_df <- add_context_seq(junc_df, size = 400, bsg = bsg)

  expect_true(nrow(cts_df) == nrow(junc_df))
  expect_true(all(c("cts_seq", "cts_junc_pos", "cts_id") %in% names(cts_df)))
})


test_that("add_context_seq is independed of junction combinations", {

  skip("only TRON local")

  # library(EnsDb.Hsapiens.v75)
  # library(tidyverse)
  # library(AnnotationDbi)

  # edb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
  # transcripts <- ensembldb::exonsBy(edb, by = "tx")

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  # transcript database
  # see: https://gitlab.rlp.net/tron/splice2neo/-/issues/38#note_235974
  # txdb_file <- "/path/to/file.txdb.sqlite"

  txdb <- loadDb(txdb_file)
  transcripts <- GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)

  # Build a GRangesList with cds composed of individual exon ranges
  # cds <- GenomicFeatures::cdsBy(txdb, by = c("tx"), use.name = TRUE)


  junc_df <- tibble::tibble(
    junc_id = c(
      "chr10_32832227_32832228_+",
      "chr10_32832227_32832228_+",
      "chr10_32832227_32832228_+"
    ),
    tx_id = c(
      "ENST00000639629.1_3",
      "ENST00000639453.1_1",
      "ENST00000639041.1_4"
    )
  )

  junc_df <- junc_df %>%
    mutate(
      tx_lst = as.list(transcripts[tx_id])
    )

  #=============================================================================
  # jx <- junc_df$junc_id %>% junc_to_gr()
  # tx_gr <- junc_df$tx_lst %>%
  #   GenomicRanges::GRangesList() %>%
  #   base::range() %>%
  #   unlist()
  # findOverlaps(jx, tx_gr)
  #=============================================================================
  # tx <- junc_df$tx_lst %>%
  #   GenomicRanges::GRangesList()
  # tx_alt <- modify_tx(tx, jx)
  #
  # SplicingGraphs::plotTranscripts(c(tx, tx_alt), from = 32750000, to = 32850000)
  # SplicingGraphs::plotTranscripts(c(tx, tx_alt))
  #=============================================================================

  cts_df <- add_context_seq(junc_df, size = 400, bsg = bsg)

  expect_true(is.na(cts_df$junc_pos_tx[[3]]))
  expect_true(is.na(cts_df$cts_seq[[3]]))

})
