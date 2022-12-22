
#' Annotate splice junctions with resulting CDS and peptide sequence
#'
#' @param df A data.frame with splice junctions in rows and at least the columns:
#'
#'   -  `junc_id` junction id consisting of genomic coordinates
#'   -  `tx_id` the ID of the affected transcript/CDS (see \code{\link{add_tx}})
#'
#' @param cds  as a named \code{\link[GenomicRanges]{GRangesList}} of coding sequences (CDS) ranges
#' @param full_pep_seq If `full_pep_seq` is set to TRUE, the output sequence will cover 15 amino acids upstream of the junction and the complete sequence until the next stop codon.
#' @param size the total size of the output sequence (might be shorter if peptide is shorter).
#' This parameter is only relevant if `full_pep_seq` is set to FALSE.
#' @param bsg \code{\link[BSgenome]{BSgenome}} object such as
#'  \code{\link[BSgenome.Hsapiens.UCSC.hg19]{BSgenome.Hsapiens.UCSC.hg19}}
#' @param keep_ranges Should GRanges of CDS and modified CDS be
#' kept? If TRUE, the list columns `cds_lst` and `cds_mod_lst` are added to the output.
#'
#' @return A data.frame as the input with the additional column(s):
#'
#'  - `cds_mod_id` an identifier made from `tx_id` and `junc_id`
#'  - `junc_pos_cds` the junction position in the modified CDS sequence
#'  - `protein` the full protein sequence of the translated modified CDS.
#'  - `protein_junc_pos` The position of the junction in the `protein` sequence
#'  - `junc_in_orf` Indicator whether the junction is located in an open reading frame.
#'  - `peptide_context_seq_raw` the peptide sequence around the junction including stop codons.
#'  - `peptide_context` the peptide sequence around the junction truncated after stop codons.
#'  - `peptide_context_junc_pos` The junction position relative to the `peptide_context` sequence
#'
#'   If the `keep_ranges` is TRUE, the following additional columns are added to
#'   the output data.frame:
#'
#'  - `cds_lst` a list of \code{\link[GenomicRanges]{GRanges}} with
#'        the original CDS as provided in `tx_id` column and `cds` object..
#'  - `cds_mod_lst` a list of \code{\link[GenomicRanges]{GRanges}} with the
#'         modified CDS ranges.
#'
#' @examples
#' requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
#' bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#'
#' add_peptide(toy_junc_df, toy_cds, size = 30, bsg = bsg)
#'
#' @export
add_peptide <- function(df, cds, full_pep_seq = TRUE, size = NULL, bsg = NULL, keep_ranges = FALSE){


  stopifnot(is.data.frame(df))
  stopifnot("junc_id" %in% names(df))
  stopifnot("tx_id" %in% names(df))
  stopifnot(class(cds) %in% c("GRangesList", "CompressedGRangesList"))
  stopifnot(is.logical(keep_ranges) & length(keep_ranges) == 1)

  # check if all input transcript IDs are in contained in the CDS object
  # stopifnot(all(df$tx_id %in% names(cds)))

  # take only the columns junc_id and tx_id and build unique combinations
  df_sub <- df %>%
    dplyr::distinct(junc_id, tx_id) %>%
    # filter for subset with matching tx_id in input cds
    dplyr::filter(tx_id %in% names(cds))


  # get GRanges as of cds
  cds_lst <- cds[df_sub$tx_id]

  if(is.null(bsg)){
    message("INFO: Use default genome sequence from BSgenome.Hsapiens.UCSC.hg19")
    bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  }

  # get junctions as GRanges object
  jx <- junc_to_gr(df_sub$junc_id)
  # test if junction is an intron retention event
  # TODO: are non IR junctions possible that follow the rule chr:pos-pos+1:strand?
  intron_retention <- ifelse(jx@ranges@width == 2, TRUE, FALSE)

  # modify transcripts by applying the splice junctions
  cds_mod <- modify_tx(cds_lst, jx)

  # get junction position in altered CDS
  junc_pos_cds = get_junc_pos(cds_mod, jx)
  # only relevant for intron retention events: get the other position of retained
  # intervall
  junc_end_tx <- get_intronretention_alt_pos(cds_lst, cds_mod, jx, intron_retention)

  # get CDS sequence
  cds_seq <- GenomicFeatures::extractTranscriptSeqs(bsg, cds_mod)

  # translate to protein seq
  suppressWarnings(
    protein <- Biostrings::translate(cds_seq, if.fuzzy.codon = "solve")
  )

  # Get context peptides around junction
  protein_junc_pos <- ceiling(junc_pos_cds / 3)
  protein_end_pos <- ceiling(junc_end_tx / 3)

  protein_len <- BiocGenerics::width(protein)


  # test if junction position is in ORF (i.e. no stop codon `*` in whole seq before)
  junc_in_orf <- XVector::subseq(protein, start = 1, end = ifelse(intron_retention & protein_end_pos < protein_junc_pos, pmin(protein_end_pos, protein_len) ,pmin(protein_junc_pos, protein_len))) %>% as.data.frame()
  junc_in_orf <- stringr::str_detect(as.character(junc_in_orf$x), "\\*", negate = TRUE)

  # extract context sequence from full peptide and cut before stop codon (*)
  if(full_pep_seq){
    pep_start <- ifelse(intron_retention & protein_end_pos < protein_junc_pos, pmax(protein_end_pos - 15 + 1, 1) ,pmax(protein_junc_pos - 15 + 1, 1))
    peptide_context_seq_raw <- XVector::subseq(protein, start = pep_start)
  }else{
    pep_start <- ifelse(intron_retention & protein_end_pos < protein_junc_pos, pmax(protein_end_pos - (size/2) + 1, 1) , pmax(protein_junc_pos - (size/2) + 1, 1))
    #pep_start <- pmax(protein_junc_pos - (size/2) + 1, 1)
    pep_end <- ifelse(intron_retention & protein_end_pos < protein_junc_pos, pmin(protein_end_pos + (size/2), protein_len) , pmin(protein_junc_pos + (size/2), protein_len))
    #pep_end <- pmin(protein_junc_pos + (size/2), protein_len)
    peptide_context_seq_raw <- XVector::subseq(protein, start = pep_start, end = pep_end)
  }

  # calculate junction position relative to context sequence
  peptide_context_junc_pos <- ifelse(intron_retention & protein_end_pos < protein_junc_pos, protein_end_pos - pep_start, protein_junc_pos - pep_start)
  #peptide_context_junc_pos <- protein_junc_pos - pep_start

  # get sequence of non-stop-codon after junction position
  peptide_context_seq_raw <- data.frame(peptide_context_seq_raw)$x
  peptide_context <- seq_truncate_nonstop(peptide_context_seq_raw, peptide_context_junc_pos)


  # Annotate table
  df_sub <- df_sub %>%
    dplyr::mutate(
      cds_mod_id = stringr::str_c(tx_id, "|", junc_id),
      junc_pos_cds = junc_pos_cds,
      protein = as.character(protein),
      protein_junc_pos = protein_junc_pos,
      junc_in_orf = junc_in_orf,

      # add NA for context sequences if the junction positon is not in an open reading frame
      peptide_context_seq_raw = ifelse(junc_in_orf, as.character(peptide_context_seq_raw), NA),
      peptide_context = ifelse(junc_in_orf, as.character(peptide_context), NA),
      peptide_context_junc_pos = ifelse(junc_in_orf, peptide_context_junc_pos, NA),
    )

  # if keep_ranges argument is TRUE add list columns of GRanges as transcripts
  if(keep_ranges){
    df_sub <- df_sub %>%
      dplyr::mutate(
        cds_lst = as.list(cds_lst),
        cds_mod_lst = as.list(cds_mod),
      )
  }

  # add annotations to input data.frame
  df <- df %>%
    dplyr::left_join(df_sub, by = c("junc_id", "tx_id"))

  return(df)

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
  suffix <- str_sub(seq, start = pos + 1)

  suffix_before_stop <- suffix %>%
    str_extract("[^*]+")

  str_c(prefix, suffix_before_stop)

}



