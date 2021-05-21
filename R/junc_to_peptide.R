
#' Annotate splice junctions with resulting CDS and peptide sequence
#'
#' @param junc_id junction id consisting of genomic coordinates
#' @param cds \code{\link[GenomicRanges]{GRangesList}} of CDS
#' @param size the total size of the output sequence (might be shorter if peptide is shorter)
#' @param bsg \code{\link[BSgenome]{BSgenome}} object such as
#'  \code{\link[BSgenome.Hsapiens.UCSC.hg19]{BSgenome.Hsapiens.UCSC.hg19}}
#'
#' @return A data.frame with annotated input junctions
#'
#' @examples
#'
#' requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
#' bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#'
#' junc_to_peptide(toy_junc_id, toy_cds, size = 30, bsg = bsg)
#'
#'
#' @export
junc_to_peptide <- function(junc_id, cds, size = 30, bsg){

  # junc_id <- toy_junc_id
  # cds <- toy_cds
  # size = 30

  #-----------------------------------------------------------------------------
  # build junction dataset with separte columns
  #-----------------------------------------------------------------------------
  junc_df <- tibble::tibble(junc_id = junc_id) %>%
    separate(junc_id, sep = "_",
             into = c("chr", "pos1", "pos2", "strand"), remove = FALSE) %>%
    mutate(
      pos1 = as.integer(pos1),
      pos2 = as.integer(pos2)
    )

  # compute affected cds
  cds_df <- junc_df %>%
    # filter out NA junctions
    filter(!is.na(junc_id)) %>%
    mutate(
      cds_df = furrr::future_pmap(
        list(chr, pos1, pos2),
        junc_to_tx, transcripts = cds)
    )

  # build junction and cds specific data sets
  cont_df <- cds_df %>%
    unnest(cds_df) %>%
    mutate(
      tx_name_alt = str_c(tx_name, "|", junc_id)
    )


  # get CDS sequence
  cds_grl <- GenomicRanges::GRangesList(cont_df$tx)
  cds_seq <- GenomicFeatures::extractTranscriptSeqs(bsg, cds_grl)

  # translate to protein seq
  suppressWarnings(
    peptide <- Biostrings::translate(cds_seq)
  )

  # Get peptides around junction
  peptide_junc_pos <- ceiling(cont_df$junc_pos_tx / 3)

  peptide_size = size
  pep_len <- BiocGenerics::width(peptide)
  pep_start <- pmax(peptide_junc_pos - (peptide_size/2) + 1, 1)
  pep_end <- pmin(peptide_junc_pos + (peptide_size/2), pep_len)

  # test if junction position is in ORF (no stop codone in whole CDS before)
  junc_in_orf <- XVector::subseq(peptide, start = 1, end = pmin(peptide_junc_pos, pep_len)) %>%
    str_detect("\\*", negate = TRUE)

  # extract context sequence from full peptide and cut before stop codon (*)
  pep_context_seq_full <- XVector::subseq(peptide, start = pep_start, end = pep_end)
  # calculate junction position relative to context sequence
  peptide_context_junc_pos <- peptide_junc_pos - pep_start

  # get sequence of non-stop-codon around junction position
  pep_context_seq_df <- seq_extract_nonstop(pep_context_seq_full,
                                            peptide_context_junc_pos)

  # Annotate table
  cont_df <- cont_df %>%
    mutate(
      peptide = as.character(peptide),
      peptide_junc_pos = peptide_junc_pos,
      junc_in_orf = junc_in_orf,
      pep_context_seq_full = as.character(pep_context_seq_full),
      peptide_context = pep_context_seq_df$seq_sub,
      peptide_context_junc_pos = pep_context_seq_df$seq_sub_pos,
    )

  junc_df %>%
    left_join(cont_df,
              by = c("junc_id", "chr", "pos1", "pos2", "strand")) %>%
    select(junc_id, tx_name, peptide, peptide_junc_pos, junc_in_orf,
           pep_context_seq_full, peptide_context, peptide_context_junc_pos)
}

#' Extract sequence at position from after previous stop codon and before next.
#'
#'  @param seq Sequence
#'  @param pos position relative to sequence
#'
#'  @return a data.frame
#'
#'  @examples
#'  seq <-  "QIP*LGSNSLLFPYQLMAGSTRP*SWALGC"
#'  seq <-  c(
#'    "QIP*LGSNSLLFPYQLMAGSTRP*SWALGC",
#'    "LKMRGDTNDILSHLD*REQRVGQ*AEAASP"
#'  )
#'  pos <- c(14, 14)
#'  seq_extract_nonstop(seq, pos)
#'
seq_extract_nonstop <- function(seq, pos){
  df <- tibble(
    seq = as.character(seq),
    pos = pos
  ) %>%
    mutate(
      cds = str_locate_all(seq, "[^*]+") %>%
        map(as_tibble),
      cds = map2(cds, pos, ~filter(.x, start <= .y, end >= .y)),
      cds = map2(cds, pos, function(cds, pos){
        if(nrow(cds) > 0)
          cds
        else
          tibble(start = pos, end = pos)
      }),
      n_cds = map_int(cds, nrow)
    ) %>%
    unnest(cds) %>%
    mutate(
      seq_sub = str_sub(seq, start, end),
      seq_sub_pos = pos - start
    )

  stopifnot(nrow(df) == length(seq))

  return(select(df, seq_sub, seq_sub_pos))

}

