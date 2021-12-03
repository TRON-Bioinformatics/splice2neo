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
  # txdb <- loadDb(txdb_file)

  gtf_url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz"
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_url)

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




test_that("add_context_seq workds for example with neg issue 46", {

  skip("only TRON local")

  # library(splice2neo)
  # library(GenomicFeatures)
  # library(AnnotationDbi)
  # library(Biostrings)

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  # transcript database
  # see: https://gitlab.rlp.net/tron/splice2neo/-/issues/38#note_235974
  # txdb_file <- "/path/to/file.txdb.sqlite"

  txdb <- loadDb(txdb_file)
  transcripts <- GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)

  df = dplyr::tibble(
      junc_id = "chr1_113153625_113162266_-",
      tx_id = "ENST00000492274.1_2"
    ) %>%
    dplyr::mutate(tx_lst = as.list(transcripts[tx_id])) %>%
    add_context_seq(size = 400, bsg = bsg)

  df$tx_lst[[1]]
  df$tx_alt_lst[[1]]
  df$cts_seq
  # "AAGAAAGCCGAGTGATGGAATCACCGCGGCGGCCGTCGGAAGGACCCGCCCGGAAACGCCGACCAAGAGGGCCCCATGTCAGTTATGCGCGGGAGTGACGTCTCCTAACGCCAGCGCCCGCGCCCGCGCCCCAGGAAGTAGGGTTTGCCTTAGATATTTGAATGGTGGTACTTCCATAAGCATGGCACATCTTTTATTGAGCAAGTATCTGTAAGCCATTTGCAACCACTGATGGGAGGAACAGAGAGCAGCATTTCAGAACCAGGTTCTCCTTCGAGGAACAGAGAAAATGAAACCAGCAGACAGAATTTGTCAGGTGACTACTTTTCTAATGTGTTTTCAGAGCTGTGTAT"
  df$cts_junc_pos

  cts_1 = stringr::str_sub(df$cts_seq, 1, df$cts_junc_pos)
  cts_2 = stringr::str_sub(df$cts_seq, df$cts_junc_pos + 1)
  # cts_1 = stringr::str_sub(df$cts_seq, 1, 153)
  # cts_2 = stringr::str_sub(df$cts_seq, 153 + 1)

  cts_1 %>% str_sub(-6)
  cts_2 %>% str_sub(1, 6)

  exp_1 <- Biostrings::getSeq(bsg, "chr1", 113153625 - 5, 113153625)
  exp_2 <- Biostrings::getSeq(bsg, "chr1", 113162266, 113162266 + 5)

  exp_1
  exp_2

  ex_gr <- GRanges(c(
    "chr1:113153620-113153625:-",
    "chr1:113162266-113162271:-"
  ))

  expected_seq = "CCTTAGATATTT"

  getSeq(bsg, ex_gr)

  # wrong order
  extractTranscriptSeqs(bsg, GRangesList(ex_gr))

  # correct order
  extractTranscriptSeqs(bsg, GRangesList(ex_gr[2:1]))

  # wrong order
  extractTranscriptSeqs(bsg, GRangesList(sort(ex_gr)))

  # fixed?
  grl <- GRangesList(ex_gr)
  grl_sorted <- S4Vectors::revElements(grl, any(strand(grl) == "-"))
  extractTranscriptSeqs(bsg, grl_sorted)

  expect_equal(as.character(extractTranscriptSeqs(bsg, grl_sorted)), expected_seq)

})
