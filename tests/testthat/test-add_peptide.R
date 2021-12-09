test_that("add_peptide works on toy example data", {

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  pep_df <- add_peptide(toy_junc_df, toy_cds, size = 30, bsg = bsg)

  expect_true(nrow(pep_df) == nrow(toy_junc_df))
  new_col_names <- c("protein", "protein_junc_pos", "peptide_context", "peptide_context_junc_pos")
  expect_true(all(new_col_names %in% names(pep_df)))

  expect_equal(unique(pep_df$junc_id), unique(toy_junc_df$junc_id))
  # check that peptides have expcted size
  expect_true(all(stringr::str_length(pep_df$peptide_context) <= 30, na.rm = TRUE))
})


test_that("add_peptide works on toy example data with keep_ranges", {

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  pep_df <- add_peptide(toy_junc_df, toy_cds, size = 30, bsg = bsg,
                            keep_ranges = TRUE)


  expect_true(nrow(pep_df) == nrow(toy_junc_df))
  new_col_names <- c("protein", "protein_junc_pos", "peptide_context", "peptide_context_junc_pos")
  expect_true(all(new_col_names %in% names(pep_df)))

  expect_equal(unique(pep_df$junc_id), unique(toy_junc_df$junc_id))
  # check that peptides have expcted size
  expect_true(all(stringr::str_length(pep_df$peptide_context) <= 30, na.rm = TRUE))
  expect_true(class(pep_df$cds_lst[[1]]) == "GRanges")
  expect_true(class(pep_df$cds_mod_lst[[1]]) == "GRanges")

})

test_that("add_peptide works when tx_id is not contained in transcripts", {

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  junc_df <- toy_junc_df
  # remove one transcript ID from transcripts object
  cds <- toy_cds[-which(names(toy_cds) == "ENST00000342992")]

  pep_df <- add_peptide(junc_df, cds, size = 30, bsg = bsg)


  expect_true(nrow(pep_df) == nrow(toy_junc_df))
  new_col_names <- c("protein", "protein_junc_pos", "peptide_context", "peptide_context_junc_pos")
  expect_true(all(new_col_names %in% names(pep_df)))
  expect_true(all(is.na(pep_df[pep_df$tx_id == "ENST00000342992", "peptide_context"])))

  # remove another transcript ==================================================
  cds <- toy_cds[-which(names(toy_cds) == "ENST00000409198")]

  pep_df <- add_peptide(junc_df, cds, size = 30, bsg = bsg)


  expect_true(nrow(pep_df) == nrow(toy_junc_df))
  new_col_names <- c("protein", "protein_junc_pos", "peptide_context", "peptide_context_junc_pos")
  expect_true(all(new_col_names %in% names(pep_df)))
  expect_true(all(is.na(pep_df[pep_df$tx_id == "ENST00000409198", "peptide_context"])))

})

test_that("add_peptide not fails for junctions outside CDS", {

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19


  # build custom junctions on chr5 which DO NOT OVERLAP with the CDS
  junc_df <- dplyr::tibble(
    junc_id = c(
      "chr5_10_20_+",
      "chr5_50_70_-"
    ),
    tx_id = names(toy_cds)[1:2]
  )

  pep_df <- add_peptide(junc_df, toy_cds, size = 30, bsg = bsg)

  expect_true(nrow(pep_df) == nrow(junc_df))
  new_col_names <- c("protein", "protein_junc_pos", "peptide_context", "peptide_context_junc_pos")
  expect_true(all(new_col_names %in% names(pep_df)))

  # expect peptide_context to be NA
  expect_true(all(is.na(pep_df$peptide_context)))

})


test_that("add_peptide works for junctions outside CDS (issue #40)", {

  # on CI avoid download of annotation data
  skip_on_ci()

  # reference data -------------------------------------------------------------
  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  gtf_url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz"
  txdb <- suppressWarnings(GenomicFeatures::makeTxDbFromGFF(gtf_url))
  cds <- GenomicFeatures::cdsBy(txdb, by = c("tx"), use.name = TRUE)
  # ----------------------------------------------------------------------------

  # chr2_220501172_220501412_+_ENST00000425141.5_1

  # build custom junctions on chr5 which DO NOT OVERLAP with the CDS
  junc_df <- dplyr::tibble(
    junc_id = "chr2_220501172_220501412_+",
    tx_id = c(
      "ENST00000425141.5_1"
    )
  )

  pep_df <- add_peptide(junc_df, cds, size = 30, bsg = bsg)

  expect_true(is.na(pep_df$peptide_context))
  expect_true(is.na(pep_df$peptide_context_junc_pos))

  # multiple junctions ---------------------------------------------------------
  # build custom junctions on chr5 which DO NOT OVERLAP with the CDS
  junc_df <- dplyr::tibble(
    junc_id = c(
      "chr2_220501172_220501411_+",
      "chr2_220501172_220501412_+",
      "chr2_220501172_220501413_+"
    ),
    tx_id = c(
      "ENST00000358055.8_4",
      "ENST00000425141.5_1",
      "ENST00000358055.8_4"
    )
  )

  pep_df <- add_peptide(junc_df, cds, size = 30, bsg = bsg)

  # expect at least one but not all peptide annotations to be NA
  expect_true(any(is.na(pep_df$peptide_context)))
  expect_true(!all(is.na(pep_df$peptide_context)))
  expect_true(any(is.na(pep_df$peptide_context_junc_pos)))
  expect_true(!all(is.na(pep_df$peptide_context_junc_pos)))

})


test_that("add_peptide works in strange combination(issue #47)", {

  # on CI avoid download of annotation data
  skip_on_ci()

  # reference data -------------------------------------------------------------
  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  gtf_url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz"
  txdb <- suppressWarnings(GenomicFeatures::makeTxDbFromGFF(gtf_url))
  cds <- GenomicFeatures::cdsBy(txdb, by = c("tx"), use.name = TRUE)
  # ----------------------------------------------------------------------------

  # use custom junctions from issue #47 See: https://gitlab.rlp.net/tron/splice2neo/-/issues/47
  # As the CDS of the first junction results in an empty range the peptide shoul be NA
  junc_enst_id <- c(
    "chr2_80526815_80531277_-|ENST00000433224.1_2", "chr2_220501172_220501412_+|ENST00000273063.10_3",
    "chr2_220501172_220501412_+|ENST00000425141.5_1", "chr2_220501172_220501412_+|ENST00000358055.8_4",
    "chr2_220501172_220501412_+|ENST00000317151.7_3"
  )
  # chr2_220501172_220501412_+_ENST00000425141.5_1

  # build custom junctions on chr5 which DO NOT OVERLAP with the CDS
  junc_df <- dplyr::tibble(
      junc_enst_id = junc_enst_id
    ) %>%
    tidyr::separate(junc_enst_id, into = c("junc_id", "tx_id"), sep = "\\|")


  pep_df <- add_peptide(junc_df, cds, size = 30, bsg = bsg)

  # As the CDS of the first junction results in an empty range the peptide shoul be NA
  # expect at least one but not all peptide annotations to be NA
  expect_true(any(is.na(pep_df$peptide_context)))
  expect_true(!all(is.na(pep_df$peptide_context)))
  expect_true(any(is.na(pep_df$peptide_context_junc_pos)))
  expect_true(!all(is.na(pep_df$peptide_context_junc_pos)))

})

test_that("seq_truncate_nonstop works on example data", {


  s1 <- seq_truncate_nonstop("1234*6789", 2) # "1234"
  expect_equal(s1, "1234")

  s2 <- seq_truncate_nonstop("1234*6789", 8) # "1234*6789"
  expect_equal(s2, "1234*6789")

  seq <- "QIP*LGSNSLLFPYQLMAGSTRP*SWALGC"
  s3 <- seq_truncate_nonstop(seq, 14) #"QIP*LGSNSLLFPYQLMAGSTRP"
  expect_equal(s3, "QIP*LGSNSLLFPYQLMAGSTRP")

})
