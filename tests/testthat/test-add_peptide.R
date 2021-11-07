test_that("add_peptide works on toy example data", {

  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
  bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

  junc_df <- tibble::tibble(
    junc_id = toy_junc_id[c(1, 6, 10)]
  ) %>%
    add_tx(toy_cds) %>%
    rename(
      cds_lst = tx_lst
    )

  pep_df <- add_peptide(junc_df, size = 30, bsg = bsg)

  expect_true(nrow(pep_df) == nrow(junc_df))
  expect_true(all(c("protein", "protein_junc_pos", "peptide_context", "peptide_context_junc_pos") %in% names(pep_df)))
})
