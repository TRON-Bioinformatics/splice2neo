test_that("add_peptide works on toy example data", {

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  junc_df <- toy_junc_df %>%
    mutate(
      cds_lst = as.list(toy_cds[tx_id])
    )

  pep_df <- add_peptide(junc_df, size = 30, bsg = bsg)

  expect_true(nrow(pep_df) == nrow(junc_df))
  new_col_names <- c("protein", "protein_junc_pos", "peptide_context", "peptide_context_junc_pos")
  expect_true(all(new_col_names %in% names(pep_df)))

  expect_equal(unique(pep_df$junc_id), unique(junc_df$junc_id))
  # check that peptides have expcted size
  expect_true(all(stringr::str_length(pep_df$peptide_context) <= 30, na.rm = TRUE))
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
