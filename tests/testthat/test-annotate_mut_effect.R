test_that("annotate_mut_effect works on toy example", {

  spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")
  df_raw <- parse_spliceai(spliceai_file)
  df <- format_spliceai(df_raw)

  annot_df <- annotate_mut_effect(df, toy_transcripts, toy_transcripts_gr)

  expect_true(nrow(annot_df) >= nrow(df))
  expect_true(length(unique(annot_df$tx_id)) > 1)

})


test_that("annotate_mut_effect works on toy example with pangolin", {

  pangolin_file <- system.file("extdata", "spliceai_output.pangolin.vcf", package = "splice2neo")

  effect_df <- parse_pangolin(pangolin_file) %>%
    format_pangolin()

  annot_df <- annotate_mut_effect(effect_df, toy_transcripts, toy_transcripts_gr)

  expect_true(nrow(annot_df) >= nrow(effect_df))
  expect_true(length(unique(annot_df$tx_id)) > 1)

})

test_that("annotate_mut_effect works on empty tibble", {

  spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")
  df_raw <- parse_spliceai(spliceai_file)
  df_raw_empty <- df_raw %>% filter(!is.na(ID))
  df_empty <- format_spliceai(df_raw_empty)
  df_raw <- parse_spliceai(spliceai_file)
  df <- format_spliceai(df_raw)

  annot_df_empty <- annotate_mut_effect(df_empty, toy_transcripts, toy_transcripts_gr)
  annot_df <- annotate_mut_effect(df, toy_transcripts, toy_transcripts_gr)

  expect_true(nrow(annot_df_empty) == 0)
  expect_true(ncol(annot_df_empty) == ncol(annot_df))
  expect_true(all(names(annot_df_empty) == names(annot_df)))

})

test_that("annotate_mut_effect does not predict events outside of exon range", {

  # event at the end of gene/transcript with only one exon
  # test that no prediction of an "intron retention"
  df <- dplyr::tibble(
    mut_id = "chr2_152043101_G_T",
    mut_effect_id = str_c(mut_id, "_", 1),
    effect = "DL",
    prob = 0.32,
    pos_rel = 8,
    chr ="chr2",
    pos = 152043101 + 8
  )

  annot_df <- annotate_mut_effect(df, toy_transcripts, toy_transcripts_gr)

  expect_true(nrow(annot_df) == 0 )


})


test_that("annotate_mut_effect works for multiple effects from same mutation", {

  skip("Long download")

  gtf_url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz"

  # parse GTF file as txdb object
  txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_url)
  transcripts <- GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)
  transcripts_gr <- GenomicFeatures::transcripts(txdb)


  df <- dplyr::tibble(
    mut_id = c(
      "chr2_62934431_G_T",
      "chr2_62934431_G_T"
    ),
    chr = c("chr2", "chr2"),
    pos_rel = c(16L, -1L),
    pos = c(62934431, 62934431) + pos_rel,
    effect = structure(3:4, .Label = c("AG","AL", "DG", "DL"), class = "factor"),
    prob = c(0.32, 0.95)
  )

  annot_df <- annotate_mut_effect(df, transcripts, transcripts_gr)

  # expect_true(nrow(annot_df) >= nrow(df))

})


test_that("annotate_mut_effect works for donor gain and aceptor gain in introns", {

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
    ~pos, ~effect, ~name,
    7,    "DG",  "DG exon",
    17,   "AG",  "AG exon",
    10,   "DG",  "DG intron",
    12,   "AG",  "AG intron",
    23,   "DL",  "DL",
    15,   "AL",  "AL",
    19,   "DL",  "DL exon"
    ) %>%
    mutate(chr = "1", pos_rel = 0, REF = "G", ALT = "T")

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

  annot_df <- annotate_mut_effect(toy_df, toy_tx, range(toy_tx))

  expect_equal(sort(annot_df$junc_id), sort(expected_jx))

})


test_that("annotate_mut_effect does not annotate intron-retention at positions within introns", {

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
    ~pos, ~effect, ~name,
    10,   "DL",  "DL intron",
    12,   "AL",  "AL intron",
    10,   "DG",  "DG intron",
    12,   "AG",  "AG intron",
    23,   "DL",  "DL",
    15,   "AL",  "AL",
  ) %>%
    mutate(chr = "1", pos_rel = 0, REF = "G", ALT = "T")

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

  annot_df <- annotate_mut_effect(toy_df, toy_tx, range(toy_tx))

  expect_equal(sort(annot_df$junc_id), sort(expected_jx))

})

