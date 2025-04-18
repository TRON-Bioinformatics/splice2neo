test_that("add_peptide works on toy example data", {

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  pep_df <- add_peptide(toy_junc_df, toy_cds, bsg = bsg)

  expect_true(nrow(pep_df) == nrow(toy_junc_df))
  new_col_names <- c("protein", "protein_junc_pos", "peptide_context", "peptide_context_junc_pos")
  expect_true(all(new_col_names %in% names(pep_df)))

  expect_equal(unique(pep_df$junc_id), unique(toy_junc_df$junc_id))
  # check that peptides have expcted size
  expect_true(all(stringr::str_length(pep_df$peptide_context) <= 35, na.rm = TRUE))
})

test_that("flanking_size parameter in add_peptide parameter works ", {

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  pep_df <- add_peptide(toy_junc_df, toy_cds, bsg = bsg, flanking_size = 14)

  test_junc <- pep_df %>% filter(junc_id == "chr2:152389996-152392205:-")
  expect_true(nchar(test_junc$peptide_context) == 16)
  expect_true(test_junc$peptide_context_junc_pos == 14)

  pep_df <- add_peptide(toy_junc_df, toy_cds, bsg = bsg, flanking_size = 13)

  test_junc <- pep_df %>% filter(junc_id == "chr2:152389996-152392205:-")
  expect_true(nchar(test_junc$peptide_context) == 15)
  expect_true(test_junc$peptide_context_junc_pos == 13)

  test_junc <- pep_df %>% filter(junc_id == "chr2:152388410-152392205:-")
  expect_true(nchar(test_junc$peptide_context) == 26)
  expect_true(test_junc$peptide_context_junc_pos == 13)
  # left side
  expect_true(substr(test_junc$peptide_context, 1, test_junc$peptide_context_junc_pos) == "NRHFKYATQLMNE")
  # right side
  expect_true(substr(test_junc$peptide_context, test_junc$peptide_context_junc_pos + 1, nchar(test_junc$peptide_context) ) == "IKYRKNYEKSKDK")

})


test_that("add_peptide_seq does not fail on predicted intron retentions at the end and beginning of a transcript ", {

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  df <- tibble(junc_id = c("chr2:152236047-152236048", "chr2:152214180-152214181"),
               tx_id = c("ENST00000243347","ENST00000243347" ))

  wt_pep_seq <- "MIILIYLFLLLWEDTQGWGFKDGIFHNSIWLERAAGVYHREARSGKYKLTYAEAKAVCEFEGGHLATYKQLEAARKIGFHVCAAGWMAKGRVGYPIVKPGPNCGFGKTGIIDYGIRLNRSERWDAYCYNPHAKECGGVFTDPKQIFKSPGFPNEYEDNQICYWHIRLKYGQRIHLSFLDFDLEDDPGCLADYVEIYDSYDDVHGFVGRYCGDELPDDIISTGNVMTLKFLSDASVTAGGFQIKYVAMDPVSKSSQGKNTSTTSTGNKNFLAGRFSHL"

  pep_df <- add_peptide(df, toy_cds, bsg = bsg)

  expect_true(nrow(pep_df) == nrow(df))
  expect_true(all(is.na(pep_df$peptide_context_seq_raw)))
  expect_true(gsub("\\*","",pep_df$protein[1]) == wt_pep_seq)

})

test_that("add_peptide works for IRs", {

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  df <- dplyr::tibble(
    junc_id = c(
      "chr17:41267742-41267743:-", # junc not in orf
      "chr2:166451723-166451724:+",
      "chr2:166533118-166533119:+", # IR
      "chr2:166535210-166535211:+", #this defines same IR as previous junc
      "chr2:166165175-166165176:+"# junc not in orf
    ),
    tx_id = c("ENST00000642945", "ENST00000314499", "ENST00000314499", "ENST00000314499", "ENST00000486878")
  )


  pep_df <- add_peptide(df, toy_cds, bsg = bsg)
  # "chr2:166535210-166535211:+" and "chr2:166535210-166535211:+" two juncs of same IR --> should be same pep seq
  expect_true(pep_df$peptide_context[3] == pep_df$peptide_context[4])
  expect_true(all(!is.na(pep_df$protein)))
  expect_true(all(!is.na(pep_df$protein)))

  pep_na <- pep_df %>% filter(str_detect(pattern = protein, string = protein_wt )) %>% pull(peptide_context)
  expect_true(all(is.na(pep_na)))

})

test_that("add_peptide returns expected results", {

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  df <- dplyr::tibble(
    junc_id = c(
      "chr2:166451723-166514269:+", # AS event
      "chr2:166451726-166514271:+", # AS event in frame
      "chr2:166451725-166514271:+", # AS event out frame
      "chr2:166451723-166532822:+", # exon skipping event
      "chr2:166533118-166533119:+" # IR event
    ),
    tx_id = c( "ENST00000314499", "ENST00000314499", "ENST00000314499", "ENST00000314499", "ENST00000314499")
  )


  pep_df <- add_peptide(df, toy_cds, bsg = bsg)

  expect_true(pep_df$peptide_context[1] == "SGDSVNPSTSSHFTQLPPFSKGRND")
  expect_true(pep_df$peptide_context[3] == "SGDSVNPSTSSHFTRLPPFSKGRND")
  expect_true(nchar(pep_df$peptide_context[2]) > nchar(pep_df$peptide_context[3]))
  expect_true(pep_df$peptide_context[5] == "PDTCTCSLAGIKCQVRVGNSGHPKTRHCLPPPKEAVMVPIMKLPTCKRKFFGKARHVIQEERHNSGREKD")

})

test_that("add_peptide works on toy example data with keep_ranges", {

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  pep_df <- add_peptide(toy_junc_df, toy_cds, bsg = bsg,
                            keep_ranges = TRUE)


  expect_true(nrow(pep_df) == nrow(toy_junc_df))
  new_col_names <- c("protein", "protein_junc_pos", "peptide_context", "peptide_context_junc_pos")
  expect_true(all(new_col_names %in% names(pep_df)))

  expect_equal(unique(pep_df$junc_id), unique(toy_junc_df$junc_id))
  # check that peptides have expcted size
  expect_true(all(stringr::str_length(pep_df$peptide_context) > 14, na.rm = TRUE))
  expect_true(class(pep_df$cds_lst[[1]]) == "GRanges")
  expect_true(class(pep_df$cds_mod_lst[[1]]) == "GRanges")

})

test_that("add_peptide works when tx_id is not contained in transcripts", {

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  junc_df <- toy_junc_df
  # remove one transcript ID from transcripts object
  cds <- toy_cds[-which(names(toy_cds) == "ENST00000342992")]

  pep_df <- add_peptide(junc_df, cds,  bsg = bsg)


  expect_true(nrow(pep_df) == nrow(toy_junc_df))
  new_col_names <- c("protein", "protein_junc_pos", "peptide_context", "peptide_context_junc_pos")
  expect_true(all(new_col_names %in% names(pep_df)))
  expect_true(all(is.na(pep_df[pep_df$tx_id == "ENST00000342992", "peptide_context"])))

  # remove another transcript ==================================================
  cds <- toy_cds[-which(names(toy_cds) == "ENST00000409198")]

  pep_df <- add_peptide(junc_df, cds, bsg = bsg)


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
      "chr5:10-20:+",
      "chr5:50-70:-"
    ),
    tx_id = names(toy_cds)[1:2]
  )

  pep_df <- add_peptide(junc_df, toy_cds, bsg = bsg)

  expect_true(nrow(pep_df) == nrow(junc_df))
  new_col_names <- c("protein", "protein_junc_pos", "peptide_context", "peptide_context_junc_pos")
  expect_true(all(new_col_names %in% names(pep_df)))

  # expect peptide_context to be NA
  expect_true(all(is.na(pep_df$peptide_context)))

})


test_that("add_peptide works for junctions outside CDS (issue #40)", {

  # on CI avoid download of annotation data
  skip("Long download")

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
    junc_id = "chr2:220501172-220501412:+",
    tx_id = c(
      "ENST00000425141.5_1"
    )
  )

  pep_df <- add_peptide(junc_df, cds,  bsg = bsg)

  expect_true(is.na(pep_df$peptide_context))
  expect_true(is.na(pep_df$peptide_context_junc_pos))

  # multiple junctions ---------------------------------------------------------
  # build custom junctions on chr5 which DO NOT OVERLAP with the CDS
  junc_df <- dplyr::tibble(
    junc_id = c(
      "chr2:220501172-220501411:+",
      "chr2:220501172-220501412:+",
      "chr2:220501172-220501413:+"
    ),
    tx_id = c(
      "ENST00000358055.8_4",
      "ENST00000425141.5_1",
      "ENST00000358055.8_4"
    )
  )

  pep_df <- add_peptide(junc_df, cds, bsg = bsg)

  # expect at least one but not all peptide annotations to be NA
  expect_true(any(is.na(pep_df$peptide_context)))
  expect_true(!all(is.na(pep_df$peptide_context)))
  expect_true(any(is.na(pep_df$peptide_context_junc_pos)))
  expect_true(!all(is.na(pep_df$peptide_context_junc_pos)))

})


test_that("add_peptide works in strange combination(issue #47)", {

  # on CI avoid download of annotation data
  skip("Long download")

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
    "chr2:80526815-80531277:-|ENST00000433224.1_2", "chr2:220501172-220501412:+|ENST00000273063.10_3",
    "chr2:220501172-220501412:+|ENST00000425141.5_1", "chr2:220501172-220501412:+|ENST00000358055.8_4",
    "chr2:220501172-220501412:+|ENST00000317151.7_3"
  )
  # chr2_220501172_220501412_+_ENST00000425141.5_1

  # build custom junctions on chr5 which DO NOT OVERLAP with the CDS
  junc_df <- dplyr::tibble(
      junc_enst_id = junc_enst_id
    ) %>%
    tidyr::separate(junc_enst_id, into = c("junc_id", "tx_id"), sep = "\\|")


  pep_df <- add_peptide(junc_df, cds, bsg = bsg)

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


test_that("add_peptide does tranlate CDS with removed start codon", {

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  # tx: ENST00000243347
  # CDS start:  chr2 152214181
  # exon1:      chr2 152214106-152214274
  # exon2:      chr2 152220457-152220594
  # junction    chr2:152214150-152220457 # removing first exon of CDS

  # exon:   ====================-------------------================
  # CDS:    ===========#########-------------------################
  # junc:        |---------------------------------|
  rmcds_junc_df <- dplyr::tibble(
    junc_id = "chr2:152214150-152220457:+",
    tx_id = "ENST00000243347"
  )

  pep_df <- add_peptide(rmcds_junc_df, toy_cds["ENST00000243347"], bsg = bsg)

  expect_false(pep_df$junc_in_orf)
  expect_true(is.na(pep_df$peptide_context))

})

test_that("is_first_reading_frame works on example data", {


  df <- dplyr::tibble(
    normalized_cds_junc_pos = c(6,7,8,9)
  )

  df1 <- df %>%
    is_first_reading_frame()

  expect_true(all(df1$is_first_reading_frame == c(TRUE, FALSE, FALSE, TRUE)))

})


test_that("annotate_junc_in_orf works on example data", {

  protein_junc_not_in_orf = "XX*XXXXXXXXXXXXXXXX"
  protein_junc_in_orf = "XXXXXXXXXXXXXX*XXXX"
  protein_junc_in_orf2 = "XXXXXXXXXXXXXXXXXX*"

  df <- dplyr::tibble(
    normalized_protein_junc_pos = c(5, 5, 5),
    protein_len = c(
      nchar(protein_junc_not_in_orf),
      nchar(protein_junc_in_orf),
      nchar(protein_junc_in_orf2)
    ),
    protein = c(protein_junc_not_in_orf, protein_junc_in_orf, protein_junc_in_orf2)
  )

  df1 <- df %>%
    annotate_junc_in_orf()

  expect_true(all(df1$junc_in_orf == c(FALSE, TRUE, TRUE)))

})


test_that("get_normalized_protein_junc_pos works for first frame and insertion of novel sequence and protein_junc_pos already on left side", {

  # toy example 1
  protein = "AAAABXXCAAAA"
  protein_wt = "AAAABCAAAA"
  junc_pos_cds = 15
  junc_pos_cds_wt = 15
  cds_length_difference = 6
  frame_shift = FALSE
  is_intron_retention = FALSE
  protein_junc_pos = 5

  df <- dplyr::tibble(
    frame_shift,
    is_intron_retention,
    protein_length_difference = nchar(protein) - nchar(protein_wt),
    cds_length_difference,
    protein,
    protein_junc_pos,
    protein_wt,
    junc_pos_cds,
    junc_pos_cds_wt
  )

  df1 <- df %>%
    get_normalized_protein_junc_pos()

  expect_true(df1$normalized_cds_junc_pos == df1$junc_pos_cds )
  expect_true(df1$normalized_protein_junc_pos == df1$protein_junc_pos )

  # toy example 2
  protein = "AAAABXXBAAAA"
  protein_wt = "AAAABCAAAA"
  junc_pos_cds = 15
  junc_pos_cds_wt = 15
  cds_length_difference = 6
  frame_shift = FALSE
  is_intron_retention = FALSE
  protein_junc_pos = 5

  df <- dplyr::tibble(
    frame_shift,
    is_intron_retention,
    protein_length_difference = nchar(protein) - nchar(protein_wt),
    protein,
    protein_junc_pos,
    cds_length_difference,
    protein_wt,
    junc_pos_cds,
    junc_pos_cds_wt
  )

  df1 <- df %>%
    get_normalized_protein_junc_pos()

  expect_true(df1$normalized_cds_junc_pos == df1$junc_pos_cds )
  expect_true(df1$normalized_protein_junc_pos == df1$protein_junc_pos )

})

test_that("get_normalized_protein_junc_pos works for first frame and insertion of novel sequence and protein_junc_pos on right side", {

  # toy example
  protein = "AAAABXXCAAAA"
  protein_wt = "AAAABCAAAA"
  junc_pos_cds = nchar("AAAABXX") * 3
  junc_pos_cds_wt = 0
  cds_length_difference = (nchar(protein) - nchar(protein_wt)) * 3
  frame_shift = FALSE
  is_intron_retention = FALSE
  protein_junc_pos = ceiling(junc_pos_cds / 3)

  df <- dplyr::tibble(
    frame_shift,
    is_intron_retention,
    cds_length_difference,
    protein,
    protein_junc_pos,
    protein_wt,
    junc_pos_cds,
    junc_pos_cds_wt
  )

  df1 <- df %>%
    get_normalized_protein_junc_pos()

  expect_true(df1$normalized_cds_junc_pos != df1$junc_pos_cds )
  expect_true(df1$normalized_protein_junc_pos != df1$protein_junc_pos )

  # TODO: check the new normalized pos!
  expect_true(df1$normalized_protein_junc_pos == 5 )

})


test_that("get_normalized_protein_junc_pos works for first frame and deletion of sequence", {

  # toy example
  protein = "AAAAAAAA"
  protein_wt = "AAAABCAAAA"
  junc_pos_cds = nchar("AAAAA") * 3
  junc_pos_cds_wt = nchar("AAAAA") * 3
  cds_length_difference = (nchar(protein) - nchar(protein_wt)) * 3
  frame_shift = FALSE
  is_intron_retention = FALSE
  protein_junc_pos = ceiling(junc_pos_cds / 3)

  df <- dplyr::tibble(
    frame_shift,
    is_intron_retention,
    protein,
    protein_junc_pos,
    cds_length_difference,
    protein_wt,
    junc_pos_cds,
    junc_pos_cds_wt
  )


  df1 <- df %>%
    get_normalized_protein_junc_pos()

  expect_true(df1$normalized_cds_junc_pos == df1$junc_pos_cds )
  expect_true(df1$normalized_protein_junc_pos == df1$protein_junc_pos )

})

test_that("get_normalized_protein_junc_pos works for non-first frame and deletion of sequence", {

  # toy example-1
  protein = "AAAAAAA"
  protein_wt = "AAAABCDAAA"
  junc_pos_cds = nchar("AAAAA") * 3 - 1
  junc_pos_cds_wt = nchar("AAAAA") * 3 - 1
  cds_length_difference = (nchar(protein) - nchar(protein_wt)) * 3
  frame_shift = FALSE
  is_intron_retention = FALSE
  protein_junc_pos = ceiling(junc_pos_cds / 3)

  df <- dplyr::tibble(
    frame_shift,
    is_intron_retention,
    protein,
    protein_junc_pos,
    cds_length_difference,
    protein_wt,
    junc_pos_cds,
    junc_pos_cds_wt
  )


  df1 <- df %>%
    get_normalized_protein_junc_pos()

  expect_true(df1$normalized_cds_junc_pos == df1$junc_pos_cds )
  expect_true(df1$normalized_protein_junc_pos != df1$protein_junc_pos )

  # toy example-2: novel AA by codon which is affected by junction
  protein = "AAADAAA"
  protein_wt = "AAAABCDAAA"
  junc_pos_cds = nchar("AAAAD") * 3 - 1
  junc_pos_cds_wt = nchar("AAAAA") * 3 - 1
  cds_length_difference = (nchar(protein) - nchar(protein_wt)) * 3
  frame_shift = FALSE
  is_intron_retention = FALSE
  protein_junc_pos = ceiling(junc_pos_cds / 3)

  df <- dplyr::tibble(
    frame_shift,
    is_intron_retention,
    protein,
    protein_junc_pos,
    cds_length_difference,
    protein_wt,
    junc_pos_cds,
    junc_pos_cds_wt
  )

  df1 <- df %>%
    get_normalized_protein_junc_pos()

  expect_true(df1$normalized_cds_junc_pos == df1$junc_pos_cds )
  expect_true(df1$normalized_protein_junc_pos != df1$protein_junc_pos )
  expect_true(df1$normalized_cds_junc_pos  ==  df1$junc_pos_cds )
  expect_true(df1$normalized_protein_junc_pos == 4 )
  expect_true(df1$exon_end_AA != df1$exon_end_AA_WT )


})

test_that("get_peptide_context works for toy data", {

  # toy example
  protein = c("AAAAAAA","AAADAAA","AAAAAAAA","AAAABXXCAAAA", "AAAABXXCAAAA", "AAAABXXBAAAA", "AAAABXXXX")
  protein_wt = c("AAAABCDAAA","AAAABCDAAA","AAAABCDAAA","AAAABCDAAA","AAAABCDAAA","AAAABCDAAA", "AAAABCDAAA")
  junc_pos_cds = c(11, 11, 12, 24, 15, 15, 15)
  junc_pos_cds_wt = c(11, 11, 12, 0, 15, 15, 15)
  cds_length_difference = c(-6, -6,-6, 6, 6, 6, -3)
  frame_shift = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE)
  is_intron_retention = FALSE
  protein_junc_pos = c(4, 4, 4, 7, 5, 5, 5)

  df <- dplyr::tibble(
    frame_shift,
    is_intron_retention,
    protein,
    protein_junc_pos,
    protein_wt,
    junc_pos_cds,
    junc_pos_cds_wt,
    protein_len = as.numeric(nchar(protein)),
    cds_length_difference,
  )

  df1 <- df %>%
    get_normalized_protein_junc_pos()

  flanking_size = 2

  df1 <- df1 %>%
    get_peptide_context(flanking_size = flanking_size)

  expect_true(all(df1$peptide_context_junc_pos == flanking_size ))
  expect_true(all(nchar(df1$peptide_context[c(1,3)]) == 2*flanking_size))
  expect_true(all(nchar(df1$peptide_context[2]) == 2*flanking_size + 1))
  expect_true(all(nchar(df1$peptide_context[4:6]) == 2*flanking_size + df1$cds_length_difference[4:6]/3 ))
  #frame-shift
  expect_true(nchar(df1$peptide_context[7]) == flanking_size + df1$protein_len[7] - df1$normalized_protein_junc_pos[7])

})

test_that("add_peptide returns NA for truncated peptides ", {

  # toy example-2
  protein = "AAAA"
  protein_wt = "AAAABCDAAA"
  junc_pos_cds = nchar("AAAA") * 3
  junc_pos_cds_wt = nchar("AAAA") * 3
  cds_length_difference = (nchar(protein) - nchar(protein_wt)) * 3
  frame_shift = FALSE
  is_intron_retention = FALSE
  protein_junc_pos = ceiling(junc_pos_cds / 3)
  protein_len = as.numeric(nchar(protein))

  df <- dplyr::tibble(
    frame_shift,
    is_intron_retention,
    protein,
    protein_junc_pos,
    cds_length_difference,
    protein_wt,
    junc_pos_cds,
    junc_pos_cds_wt,
    protein_len
  )


  df1 <- df %>%
    get_normalized_protein_junc_pos()

  df1 <- df1 %>%
    get_peptide_context(flanking_size = 2)

  expect_true(is.na(df1$peptide_context))
  expect_true(stringr::str_detect(fixed(df1$protein_wt), fixed(df1$protein)))

})




test_that("add_peptides works with custom toy example data", {

  ##############################################################################
  #       M  M  M  M  M  M  M  M  M  M
  #      <=><=><=><=><=><=><=><=><=><=>
  #seq   ATGATGATGATGATGATGATGATGATGATG
  #      0        1         2         3
  #      123456789012345678901234567890
  #cds1     ======         ======
  #jx1           |--------|
  #mcds1    ======        =======

  #cds2     ======         ======
  #jx2              |------|
  #mcds2    =========      ======

  #seq   ATGATGATGATGATGATGATGATGATGATG
  #cds3     =====          ======
  #jx3             |-------|
  #mcds3    ========       ======


  ##############################################################################

  # consturctur function to build custom geome seq.
  # Source: https://github.com/Bioconductor/BSgenome/issues/3
  build_bsgenome <- function(dna, circ_seqs=NA, genome=NA,
                       organism=NA, common_name=NA, provider=NA,
                       release_date=NA, release_name=NA, source_url=NA)
  {
    ## Some sanity checks.
    if (!is(dna, "DNAStringSet"))
      dna <- as(dna, "DNAStringSet")
    seqnames <- names(dna)
    if (is.null(seqnames))
      stop("'dna' must have names")
    if (!is.character(circ_seqs))
      circ_seqs <- as.character(circ_seqs)
    if (!identical(circ_seqs, NA_character_)) {
      if (anyNA(circ_seqs) ||
          anyDuplicated(circ_seqs) ||
          !all(circ_seqs %in% seqnames))
        stop(wmsg("when not set to NA, 'circ_seqs' must ",
                  "contain unique sequence names that are ",
                  "present in 'names(dna)'"))
    }
    stopifnot(S4Vectors::isSingleStringOrNA(genome),
              S4Vectors::isSingleStringOrNA(organism),
              S4Vectors::isSingleStringOrNA(common_name),
              S4Vectors::isSingleStringOrNA(provider),
              S4Vectors::isSingleStringOrNA(release_date),
              S4Vectors::isSingleStringOrNA(release_name),
              S4Vectors::isSingleStringOrNA(source_url))

    ## Write the sequences to disk.
    seqs_dirpath <- tempfile(pattern="BSgenome_seqs_dir")
    dir.create(seqs_dirpath)
    twobit_filepath <- file.path(seqs_dirpath, "single_sequences.2bit")
    rtracklayer::export(dna, twobit_filepath)

    ## Create the BSgenome object.
    BSgenome::BSgenome(organism=as.character(organism),
                        common_name=as.character(organism),
                        provider=as.character(provider),
                        provider_version=as.character(genome),
                        release_date=as.character(release_date),
                        release_name=as.character(release_name),
                        source_url=as.character(source_url),
                        seqnames=seqnames,
                        circ_seqs=circ_seqs,
                        mseqnames=NULL,
                        seqs_pkgname=NA_character_,
                        seqs_dirpath=seqs_dirpath)
  }
  toy_bsg <- build_bsgenome(
    c(c = stringr::str_c(rep("ATG", 10), collapse = "")),
    genome="hg00")


  cds_wt <- GenomicRanges::GRangesList(list(
    cds1 = GenomicRanges::GRanges(c(
      "c:4-9:+",
      "c:19-24:+"
    )),
    cds2 = GenomicRanges::GRanges(c(
      "c:4-9:+",
      "c:19-24:+"
    )),
    cds3 = GenomicRanges::GRanges(c(
      "c:4-8:+",
      "c:19-24:+"
    ))
  ))

  jx <- GenomicRanges::GRanges(c(
    "c:9-18",
    "c:12-19",
    "c:11-19"
  )
  )

  cds_mod <- GenomicRanges::GRangesList(list(
    cds1 = GenomicRanges::GRanges(c(
      "c:4-9:+",
      "c:18-24:+"
    )),
    cds2 = GenomicRanges::GRanges(c(
      "c:4-12:+",
      "c:19-24:+"
    )),
    cds3 = GenomicRanges::GRanges(c(
      "c:4-11:+",
      "c:19-24:+"
    ))
  ))

  junc_df <- tibble(
    junc_id = as.character(jx),
    tx_id = names(cds_wt)
  )

  pep_df <- junc_df %>%
    add_peptide(cds = cds_wt, bsg = toy_bsg)

  expect_equal(
    pep_df$peptide_context,
    c("MMDD",
      "MMMMM",
      "MMI")
  )


})


test_that("annotate_truncated_cds works on toy sample ", {

  # toy example-2
  protein = c("AAAA*XXXX", "AAAAXX*XX", "AAXX*XX", "*XXXX", "", "XXXXXXX")
  protein_wt = c("AAAABCDAAA", "AAAABCDAAA", "AAAABCDAAA", "AAAABCDAAA", NA, "AAAABCDAAA")
  junc_in_orf = c(TRUE, TRUE, TRUE, TRUE, NA, FALSE)

  df <- dplyr::tibble(
    protein,
    protein_wt,
    protein_len = as.numeric(nchar(protein)),
    junc_in_orf,
    cds_length_difference = (nchar(protein) - nchar(protein_wt)) * 3
  )

  df1 <- df %>%
    annotate_truncated_cds()

  expect_true(ncol(df) + 2 == ncol(df1) )
  expect_true(df1$truncated_cds[1])
  expect_true(!df1$truncated_cds[2])
  expect_true(all(
    df1$cds_description == c(
      "truncated wt cds",
      "mutated cds",
      "mutated cds",
      "no mutated gene product",
      "no mutated gene product",
      "not in ORF"
    )
  ))

})


