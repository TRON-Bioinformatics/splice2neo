test_that("annotate_spliceai_junction works on toy example", {

  spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")
  df_raw <- parse_spliceai(spliceai_file)
  df <- format_spliceai(df_raw)

  annot_df <- annotate_spliceai_junction(df, toy_transcripts, toy_transcripts_gr)

  expect_true(nrow(annot_df) >= nrow(df))
  expect_true(length(unique(annot_df$tx_id)) > 1)

})

test_that("annotate_spliceai_junction does not predict events outside of exon range", {

  # event at the end of gene/transcript with only one exon
  # test that no prediction of an "intron retention"
  df <- dplyr::tibble(
    CHROM = c("chr2"),
    POS = c("152043101"),
    ID = c(NA_character_),
    REF = c("G"),
    ALT = c("T"),
    QUAL = c(NA_character_),
    FILTER = c(NA_character_),
    Key = c(3L),
    ALLELE = c("T"),
    change = c("DL"),
    prob = c(0.32),
    pos_rel = c(8L)
  )

  annot_df <- annotate_spliceai_junction(df, toy_transcripts, toy_transcripts_gr)

  expect_true(nrow(annot_df) == 0 )


})


test_that("annotate_spliceai_junction works for multiple effects from same mutation", {

  skip("Long download")

  gtf_url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz"

  # parse GTF file as txdb object
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_url)
  transcripts <- GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)
  transcripts_gr <- GenomicFeatures::transcripts(txdb)


  df <- dplyr::tibble(
    CHROM = c("chr2", "chr2"),
    POS = c("62934431", "62934431"),
    ID = c(NA_character_, NA_character_),
    REF = c("G", "G"),
    ALT = c("T", "T"),
    QUAL = c(NA_character_, NA_character_),
    FILTER = c(NA_character_, NA_character_),
    Key = c(3L, 3L),
    ALLELE = c("T", "T"),
    change = structure(3:4, .Label = c("AG","AL", "DG", "DL"), class = "factor"),
    prob = c(0.32, 0.95),
    pos_rel = c(16L, -1L)
  )

  annot_df <- annotate_spliceai_junction(df, transcripts, transcripts_gr)

  # expect_true(nrow(annot_df) >= nrow(df))

})


test_that("annotate_spliceai_junction works for donor gain and aceptor gain in introns", {

  #=============================================================================
  #          1         2         3         4
  # 1234567890123456789012345678901234567890
  # -#######------#########---#########-----  tx1 (positive strand)
  #       |                                   DG exon
  #       [=======]                           DG exon junction
  #                 |                         AG exon
  #        [========]                         AG exon junction
  #          |                                DG intron
  #          [====]                           DG intron junction
  #            |                              AG intron
  #        [===]                              AG intron junction
  #                       |                   DL
  #                       []                  DL intron-retention junction
  #        [==================]               DL exon-skipping junction
  #               |                           AL
  #              []                           AL intron-retention junction
  #        [==================]               AL exon-skipping junction
  #                   |                       DL exon
  # 1234567890123456789012345678901234567890
  #          1         2         3         4
  # -#######------#########---#########-----  tx2 (negative strand)
  #       |                                   DG exon
  #                                           DG exon junction
  #                 |                         AG exon
  #                 [=========]               AG exon junction
  #          |                                DG intron
  #        [=]                                DG intron junction
  #            |                              AG intron
  #            [==]                           AG intron junction
  #                       |                   DL
  #               |                           AL
  #                   |                       DL exon
  #=============================================================================

  toy_df <- dplyr::tribble(
    ~POS, ~change, ~name,
    7,    "DG",  "DG exon",
    17,   "AG",  "AG exon",
    10,   "DG",  "DG intron",
    12,   "AG",  "AG intron",
    23,   "DL",  "DL",
    15,   "AL",  "AL",
    19,   "DL",  "DL exon"
    ) %>%
    mutate(CHROM = "1", pos_rel = 0, REF = "G", ALT = "T")

  toy_tx <- GenomicRanges::GRangesList(list(
    tx1 = GenomicRanges::GRanges(
      c("1", "1", "1"),
      IRanges::IRanges(
        c(2, 15, 27),
        c(8, 23, 35),
      ),
      strand = c("+", "+", "+")
    ),
    tx2 = GenomicRanges::GRanges(
      c("1", "1", "1"),
      IRanges::IRanges(
        c(2, 15, 27),
        c(8, 23, 35),
      ),
      strand = c("-", "-", "-")
    )
  ))

  expected_pos = c("1:7-15:+", "1:8-17:+", "1:10-15:+", "1:8-12:+", "1:23-24:+",
                  "1:8-27:+", "1:14-15:+", "1:8-27:+")
  expected_neg = c("1:17-27:-", "1:8-10:-", "1:12-15:-")
  expected_jx = c(expected_pos, expected_neg)

  annot_df <- annotate_spliceai_junction(toy_df, toy_tx, range(toy_tx))

  expect_equal(sort(annot_df$junc_id), sort(expected_jx))

})


test_that("annotate_spliceai_junction does not annotate intron-retention at positions within introns", {

  #=============================================================================
  #          1         2         3         4
  # 1234567890123456789012345678901234567890
  # -#######------#########---#########-----  tx1 (positive strand)
  #          |                                DL intron
  #            |                              AL intron
  #          |                                DG intron
  #          [====]                           DG intron junction
  #            |                              AG intron
  #        [===]                              AG intron junction
  #                       |                   DL
  #                       []                  DL intron-retention junction
  #        [==================]               DL exon-skipping junction
  #               |                           AL
  #              []                           AL intron-retention junction
  #        [==================]               AL exon-skipping junction
  # 1234567890123456789012345678901234567890
  #          1         2         3         4
  # -#######------#########---#########-----  tx2 (negative strand)
  #          |                                DL intron
  #            |                              AL intron
  #          |                                DG intron
  #        [=]                                DG intron junction
  #            |                              AG intron
  #            [==]                           AG intron junction
  #                       |                   DL
  #               |                           AL
  #=============================================================================

  toy_df <- dplyr::tribble(
    ~POS, ~change, ~name,
    10,   "DL",  "DL intron",
    12,   "AL",  "AL intron",
    10,   "DG",  "DG intron",
    12,   "AG",  "AG intron",
    23,   "DL",  "DL",
    15,   "AL",  "AL",
  ) %>%
    mutate(CHROM = "1", pos_rel = 0, REF = "G", ALT = "T")

  toy_tx <- GenomicRanges::GRangesList(list(
    tx1 = GenomicRanges::GRanges(
      c("1", "1", "1"),
      IRanges::IRanges(
        c(2, 15, 27),
        c(8, 23, 35),
      ),
      strand = c("+", "+", "+")
    ),
    tx2 = GenomicRanges::GRanges(
      c("1", "1", "1"),
      IRanges::IRanges(
        c(2, 15, 27),
        c(8, 23, 35),
      ),
      strand = c("-", "-", "-")
    )
  ))

  expected_pos = c("1:10-15:+", "1:8-12:+", "1:23-24:+",
                   "1:8-27:+", "1:14-15:+", "1:8-27:+")
  expected_neg = c("1:8-10:-", "1:12-15:-")
  expected_jx = c(expected_pos, expected_neg)

  annot_df <- annotate_spliceai_junction(toy_df, toy_tx, range(toy_tx))

  expect_equal(sort(annot_df$junc_id), sort(expected_jx))

})

