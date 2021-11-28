
#' Annotate splice junctions with resulting CDS and peptide sequence
#'
#' @param df A data.frame with splice junctions in rows and at least the columns:
#'
#'   -  `junc_id` junction id consisting of genomic coordinates
#'   -  `tx_id` the ID of the affected transcript (see \code{\link{add_tx}})
#'   -  `cds_lst` a list of \code{\link[GenomicRanges]{GRanges}} with the
#'         CDS (see \code{\link{add_tx}})
#'
#' @param size the total size of the output sequence (might be shorter if peptide is shorter)
#' @param bsg \code{\link[BSgenome]{BSgenome}} object such as
#'  \code{\link[BSgenome.Hsapiens.UCSC.hg19]{BSgenome.Hsapiens.UCSC.hg19}}
#'
#' @return A data.frame as the input with the additional column(s):
#'
#'  - `cds_alt_lst` a list of \code{\link[GenomicRanges]{GRanges}} with the
#'         modified CDS ranges.
#'  - `cds_id_alt` an identifier made from `tx_id` and `junc_id`
#'  - `junc_pos_cds` the junction position in the modified CDS sequence
#'  - `protein` the full protein sequence of the translated modified CDS.
#'  - `protein_junc_pos` The position of the junction in the `protein` sequence
#'  - `junc_in_orf` Indicator whether the junction is located in an open reading frame.
#'  - `peptide_context_seq_raw` the peptide sequence around the junction including stop codons.
#'  - `peptide_context` the peptide sequence around the junction truncated after stop codons.
#'  - `peptide_context_junc_pos` The junction position relative to the `peptide_context` sequence
#'
#' @examples
#' requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
#' bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#'
#' junc_df <- toy_junc_df %>%
#'   dplyr::mutate(
#'     cds_lst = as.list(toy_cds[tx_id])
#'   )
#'
#' add_peptide(junc_df, size = 30, bsg = bsg)
#'
#' @export
add_peptide <- function(df, size = 30, bsg = NULL){


  stopifnot(is.data.frame(df))
  stopifnot("junc_id" %in% names(df))
  stopifnot("tx_id" %in% names(df))
  stopifnot("cds_lst" %in% names(df))

  if(is.null(bsg)){
    message("INFO: Use default genome sequence from BSgenome.Hsapiens.UCSC.hg19")
    bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  }

  # get junctions as GRanges object
  jx <- junc_to_gr(df$junc_id)

  # modify transcripts by applying the splice junctions
  cds_alt <- modify_tx(df$cds_lst, jx)

  # build new id
  cds_id_alt = str_c(df$tx_id, "|", df$junc_id)

  # get junction position in altered CDS
  junc_pos_cds = get_junc_pos(cds_alt, jx)

  # get CDS sequence
  cds_seq <- GenomicFeatures::extractTranscriptSeqs(bsg, cds_alt)

  # translate to protein seq
  suppressWarnings(
    protein <- Biostrings::translate(cds_seq, if.fuzzy.codon = "solve")
  )

  # Get context peptides around junction
  protein_junc_pos <- ceiling(junc_pos_cds / 3)

  protein_len <- BiocGenerics::width(protein)
  pep_start <- pmax(protein_junc_pos - (size/2) + 1, 1)
  pep_end <- pmin(protein_junc_pos + (size/2), protein_len)


  # test if junction position is in ORF (i.e. no stop codon `*` in whole seq before)
  junc_in_orf <- XVector::subseq(protein, start = 1, end = pmin(protein_junc_pos, protein_len)) %>%
    str_detect("\\*", negate = TRUE)

  # extract context sequence from full peptide and cut before stop codon (*)
  peptide_context_seq_raw <- XVector::subseq(protein, start = pep_start, end = pep_end)

  # calculate junction position relative to context sequence
  peptide_context_junc_pos <- protein_junc_pos - pep_start

  # get sequence of non-stop-codon after junction position
  peptide_context <- seq_truncate_nonstop(peptide_context_seq_raw, peptide_context_junc_pos)


  # Annotate table
  df %>%
    dplyr::mutate(
      cds_alt_lst = as.list(cds_alt),
      cds_id_alt = cds_id_alt,
      junc_pos_cds = junc_pos_cds,
      protein = as.character(protein),
      protein_junc_pos = protein_junc_pos,
      junc_in_orf = junc_in_orf,

      # add NA for context sequences if the junction positon is not in an open reading frame
      peptide_context_seq_raw = ifelse(junc_in_orf, as.character(peptide_context_seq_raw), NA),
      peptide_context = ifelse(junc_in_orf, as.character(peptide_context), NA),
      peptide_context_junc_pos = ifelse(junc_in_orf, peptide_context_junc_pos, NA),
    )

}

#' Truncate input sequence after input position before next stop codons (`*`).
#'
#' @param seq Sequence
#' @param pos position relative to sequence
#' @return a character
#'
#' This can be used to remove peptide sequence parts in a context sequence (in
#' a fixed-sized window around a position of interest) after a stop codon `*`
#' occurs.
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



