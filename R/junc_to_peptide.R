
#' Annotate splice junctions with resulting CDS and peptide sequence
#'
#' @param junc_id junction id consisting of genomic coordinates
#' @param cds \code{\link[GenomicRanges]{GRangesList}} of CDS
#' @param tx_id a vector of transcript ids in the same length as `junc_id`. If
#'    provided, only transcripts with the associated id will be considered for
#'    each junction. If NA (default) no input CDS will be considered.
#' @param size the total size of the output sequence (might be shorter if peptide is shorter)
#' @param bsg \code{\link[BSgenome]{BSgenome}} object such as
#'  \code{\link[BSgenome.Hsapiens.UCSC.hg19]{BSgenome.Hsapiens.UCSC.hg19}}
#'
#' @return A data.frame with annotated input junctions
#'
#'  - `junc_id` the input `junc_id`
#'  - `tx_id_input` the input `tx_id`
#'  - `tx_id` the id of the used cds
#'  - `tx_id_alt` a combination of the used cds id with the junction id
#'  - `peptide` the full pepetide sequence of the translated cds.
#'  - `peptide_junc_pos` The position of the junctino in the `peptide` sequence
#'  - `junc_in_orf` Indicator whehter the junction is located in an open reading frame.
#'  - `pep_context_seq_full` the peptide sequence around the junction includig stop codons.
#'  - `peptide_context` the peptide sequence around the junction without stop codons.
#'  - `peptide_context_junc_pos` The junction position relative to the `peptide_context` sequence
#'
#' @examples
#'
#' requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
#' bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#'
#' junc_to_peptide("chr2_152389996_152392205_-", toy_cds, size = 30, bsg = bsg)
#'
#' junc_to_peptide(toy_junc_id, toy_cds, tx_id = toy_junc_id_enst, size = 30, bsg = bsg)
#'
#'@import dplyr
#'
#' @export
junc_to_peptide <- function(junc_id, cds, tx_id = NA, size = 30, bsg = NULL){

  if(is.null(bsg)){
    message("INFO: Use default genome sequence from BSgenome.Hsapiens.UCSC.hg19")
    bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  }

  # junc_id <- toy_junc_id
  # tx_id <- toy_junc_id_enst
  # cds <- toy_cds
  # size = 30

  #-----------------------------------------------------------------------------
  # build junction dataset with separte columns
  #-----------------------------------------------------------------------------
  junc_df <- tibble::tibble(
    junc_id = junc_id,
    tx_id_input = tx_id
    ) %>%
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

      # get the relevant cds, if ID given use specific cds, otherwise take all
      cds_grl = furrr::future_map_if(
        .x = tx_id_input,
        .p = ~!is.na(.x),
        .f = ~cds[.],
        .else = ~ cds
      ),

      # apply junc_to_tx() to modify transcript ranges according to the junction
      cds_df = furrr::future_pmap(
        list(chr, pos1, pos2, cds_grl),
        junc_to_tx)
    )

  # build junction and cds specific data sets
  cont_df <- cds_df %>%
    unnest(cds_df) %>%
    mutate(
      tx_id_alt = str_c(tx_id, "|", junc_id)
    )


  # get CDS sequence
  if(nrow(cont_df) > 0){
    cds_grl <- GenomicRanges::GRangesList(cont_df$tx)
  }else{
    cds_grl <- GenomicRanges::GRangesList()
  }
  cds_seq <- GenomicFeatures::extractTranscriptSeqs(bsg, cds_grl)

  # translate to protein seq
  suppressWarnings(
    peptide <- Biostrings::translate(cds_seq)
  )

  # Get peptides around junction
  if(nrow(cont_df) > 0){
    peptide_junc_pos <- ceiling(cont_df$junc_pos_tx / 3)
  }else{
    peptide_junc_pos <- numeric()
  }

  pep_len <- BiocGenerics::width(peptide)
  pep_start <- pmax(peptide_junc_pos - (size/2) + 1, 1)
  pep_end <- pmin(peptide_junc_pos + (size/2), pep_len)

  # test if junction position is in ORF (no stop codone in whole CDS before)
  junc_in_orf <- XVector::subseq(peptide, start = 1, end = pmin(peptide_junc_pos, pep_len)) %>%
    str_detect("\\*", negate = TRUE)

  # extract context sequence from full peptide and cut before stop codon (*)
  pep_context_seq_full <- XVector::subseq(peptide, start = pep_start, end = pep_end)
  # calculate junction position relative to context sequence
  peptide_context_junc_pos <- peptide_junc_pos - pep_start

  # get sequence of non-stop-codon after junction position
  peptide_context <- seq_truncate_nonstop(pep_context_seq_full, peptide_context_junc_pos)

  # Annotate table
  cont_df <- cont_df %>%
    mutate(
      peptide = as.character(peptide),
      peptide_junc_pos = peptide_junc_pos,
      junc_in_orf = junc_in_orf,
      pep_context_seq_full = as.character(pep_context_seq_full),
      peptide_context = peptide_context,
      peptide_context_junc_pos = peptide_context_junc_pos,
    )

  junc_df %>%
    left_join(cont_df,
              by = c("junc_id", "chr", "pos1", "pos2", "strand", "tx_id_input")) %>%
    select(junc_id, tx_id_input, tx_id, tx_id_alt, peptide, peptide_junc_pos, junc_in_orf,
           pep_context_seq_full, peptide_context, peptide_context_junc_pos)

}


#' Truncate input sequence after input position before next stop codons (`*`).
#'
#' @param seq Sequence
#' @param pos position relative to sequence
#' @return a character
#'
#' @examples
#'
#' seq_truncate_nonstop("1234*6789", 2) # "1234"
#' seq_truncate_nonstop("1234*6789", 8) # "1234*6789"
#'
#' seq <- "QIP*LGSNSLLFPYQLMAGSTRP*SWALGC"
#' seq_truncate_nonstop(seq, 14) #"QIP*LGSNSLLFPYQLMAGSTRP"
#'
#' @export
seq_truncate_nonstop <- function(seq, pos){

  prefix <- str_sub(seq, 1, pos)
  suffix <-  str_sub(seq, start = pos + 1)

  suffix_before_stop <- suffix %>%
    str_extract("[^*]+")

  str_c(prefix, suffix_before_stop)

}


#' Extract sub-sequence around input position between stop codons (`*`).
#'
#' @param seq Sequence
#' @param pos position relative to sequence
#' @return a data.frame
#'
#' @examples
#'
#' seq <- "QIP*LGSNSLLFPYQLMAGSTRP*SWALGC"
#' seq <- c(
#'  "QIP*LGSNSLLFPYQLMAGSTRP*SWALGC",
#'  "LKMRGDTNDILSHLD*REQRVGQ*AEAASP"
#' )
#' pos <- c(14, 14)
#' seq_extract_nonstop(seq, pos)
#'
#' @export
seq_extract_nonstop <- function(seq, pos){
  df <- tibble(
    seq = as.character(seq),
    pos = pos
  ) %>%
    mutate(
      cds = str_locate_all(seq, "[^*]+") %>%
        purrr::map(as_tibble),
      cds = purrr::map2(cds, pos, ~filter(.x, start <= .y, end >= .y)),
      cds = purrr::map2(cds, pos, function(cds, pos){
        if(nrow(cds) > 0)
          cds
        else
          tibble(start = pos, end = pos)
      }),
      n_cds = purrr::map_int(cds, nrow)
    ) %>%
    unnest(cds)

  if(nrow(df) > 0){
    df <- df %>%
      mutate(
        seq_sub = str_sub(seq, start, end),
        seq_sub_pos = pos - start
      )
  }else{
    df <- df %>%
      mutate(
        seq_sub = character(),
        seq_sub_pos = numeric()
      )

  }

  stopifnot(nrow(df) == length(seq))

  return(select(df, seq_sub, seq_sub_pos))

}



