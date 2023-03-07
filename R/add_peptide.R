
#' Annotate splice junctions with resulting CDS and peptide sequence
#'
#' @param df A data.frame with splice junctions in rows and at least the columns:
#'
#'   -  `junc_id` junction id consisting of genomic coordinates
#'   -  `tx_id` the ID of the affected transcript/CDS (see \code{\link{add_tx}})
#'
#' @param cds  as a named \code{\link[GenomicRanges]{GRangesList}} of coding sequences (CDS) ranges
#' @param flanking_size the number of wild-type amino acids flanking the junction or novel sequence caused by the splice junction to the left and to the right.
#' Frame shift junctions are flanked by wild-type amino acids to the left and are annotated until the first stop codon.
#' @param bsg \code{\link[BSgenome]{BSgenome}} object such as
#'  \code{\link[BSgenome.Hsapiens.UCSC.hg19]{BSgenome.Hsapiens.UCSC.hg19}}
#' @param keep_ranges Should GRanges of CDS and modified CDS be
#' kept? If TRUE, the list columns `cds_lst` and `cds_mod_lst` are added to the output.
#'
#' @return A data.frame as the input with the additional column(s):
#'
#'  - `cds_mod_id` an identifier made from `tx_id` and `junc_id`
#'  - `junc_pos_cds` the junction position in the modified CDS sequence
#'  - `frame_shift` Indicator whether junction leads to frame shift.
#'  - `is_first_reading_frame` Indicator whether junction is in 1st reading frame for in-frame junctions.
#'  - `normalized_cds_junc_pos` The normalized position of the junction in the modified CDS sequence to the left junction side.
#'  - `protein` the full protein sequence of the translated modified CDS.
#'  - `normalized_protein_junc_pos` The normalized position of the junction in the `protein` sequence to the left junction side.
#'  - `peptide_context_junc_pos` The junction position relative to the `peptide_context` sequence
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
#' add_peptide(toy_junc_df, toy_cds, bsg = bsg, flanking_size = 13)
#'
#' @export
add_peptide <- function(df, cds, flanking_size = 14, bsg = NULL, keep_ranges = FALSE){


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
  intron_retention <- jx@ranges@width == 2

  # modify transcripts by applying the splice junctions
  cds_mod <- modify_tx(cds_lst, jx)

  # get junction position in altered CDS
  junc_pos_cds = get_junc_pos(cds_mod, jx)
  junc_pos_cds_wt = get_junc_pos(cds_lst, jx)
  junc_in_cds = !is.na(junc_pos_cds)

  # get CDS sequence
  cds_seq <- GenomicFeatures::extractTranscriptSeqs(bsg, cds_mod)
  cds_seq_wt <- GenomicFeatures::extractTranscriptSeqs(bsg, cds_lst)

  # frame shift in modified cds seq?
  cds_length_difference  <-
    (BiocGenerics::width(cds_seq) - BiocGenerics::width(cds_seq_wt))
  frame_shift <- cds_length_difference %% 3 > 0

  # translate to protein seq
  suppressWarnings(
    protein <- Biostrings::translate(cds_seq, if.fuzzy.codon = "solve")
  )
  suppressWarnings(
    protein_wt <- Biostrings::translate(cds_seq_wt, if.fuzzy.codon = "solve")
  )

  # extract context sequence from full peptide and cut before stop codon (*)
  df_positions <- df_sub %>%
    dplyr::mutate(
    intron_retention = intron_retention,
    strand = str_sub(df_sub$junc_id, -1),
    protein = protein %>% as.character(),
    protein_wt = protein_wt %>% as.character(),
    frame_shift = frame_shift,
    junc_in_cds = junc_in_cds,
    cds_mod_id = stringr::str_c(tx_id, "|", junc_id),
    cds_length_difference = cds_length_difference,
    junc_pos_cds = junc_pos_cds,
    junc_pos_cds_wt = junc_pos_cds_wt,
    # Get context peptides around junction
    protein_junc_pos = ceiling(junc_pos_cds / 3),
    protein_junc_pos_not = junc_pos_cds / 3,
    # end for IRs
    protein_length_difference = ifelse(!frame_shift &
                                         !intron_retention, cds_length_difference / 3, NA),
    protein_len = as.numeric(BiocGenerics::width(protein)),
  )

  df_positions <- df_positions %>%
    is_first_reading_frame() %>%
    get_normalized_protein_junc_pos()%>%
    annotate_junc_in_orf()

  df_annotated_peptide <- df_positions %>%
    get_peptide_context(flanking_size = flanking_size)

  # Annotate table
  df_annotated_peptide <- df_annotated_peptide %>%
    dplyr::mutate(
      # add NA for context sequences if the junction position is not in an open reading frame
      peptide_context_seq_raw = ifelse(junc_in_orf, as.character(peptide_context_seq_raw), NA),
      peptide_context = ifelse(junc_in_orf, as.character(peptide_context), NA),
      peptide_context_junc_pos = ifelse(junc_in_orf, peptide_context_junc_pos, NA),
    )

  # if keep_ranges argument is TRUE add list columns of GRanges as transcripts
  if(keep_ranges){
    df_annotated_peptide <- df_annotated_peptide %>%
      dplyr::mutate(
        cds_lst = as.list(cds_lst),
        cds_mod_lst = as.list(cds_mod),
      )
  }

  # df_annotated_peptide <- df_annotated_peptide %>%
  #   dplyr::select(
  #     -intron_retention,
  #     -strand,
  #     #-protein_wt,
  #     -junc_in_cds,
  #     -protein_len,
  #     -pep_start,
  #     -pep_end,
  #     -exon1_end_AA,
  #     -exon1_end_AA_WT,
  #     -exon2_start_AA_WT,
  #     -cds_length_difference,
  #     -WT_protein_length_difference,
  #     -addtional_AA,
  #     -protein_length_difference,
  #     -junc_pos_cds_wt,
  #     -protein_junc_pos_not
  #   )

  # add annotations to input data.frame
  df <- df %>%
    dplyr::left_join(df_annotated_peptide, by = c("junc_id", "tx_id"))

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


#' Determine if transcript is in first reading frame, given a tibble with:
#' `junc_pos_cds`
#'
#' @param df Sequence
#'
#' Returns `TRUE` if transcript is in first reading frame and `FALSE` if transcript follows second or third reading frame.
#' occurs.
#'
is_first_reading_frame <- function(df){

  df_mod <- df %>%
   dplyr::mutate(
     is_first_reading_frame = case_when(
       round((junc_pos_cds/3) %% 1, 2) == 0 ~ TRUE,
       round((junc_pos_cds/3) %% 1, 2) != 0 ~ FALSE,
     )
   )

}


#' Tests if junction position is in ORF
# i.e. no stop codon `*` in whole seq before junction given a tibble with:
#'  `intron_retention`
#'  `normalized_protein_junc_pos`
#'  `protein_len`
#'
#' @param df tibble
#'
#'
annotate_junc_in_orf <- function(df){

  df_mod <- df %>%
    dplyr::mutate(protein_until_junction = stringr::str_sub(
      protein,
      start = 1,
      end = pmin(normalized_protein_junc_pos, protein_len)
    )) %>%
    dplyr::mutate(junc_in_orf = stringr::str_detect(protein_until_junction, "\\*", negate = TRUE))%>%
    dplyr::mutate(junc_in_orf = ifelse(is.na(junc_in_orf), FALSE, junc_in_orf))%>%
    dplyr::select(-protein_until_junction)

  return(df_mod)
}


#' Annotate the normalized junction position in the protein,
# i.e. the postion of the last WT amino acid in the mutated protein on the left side
#'  `frame_shift`
#'  `protein_length_difference`
#'  `protein_junc_pos`
#'  `protein`
#'  `protein_wt`
#'  `junc_pos_cds`
#'
#' @param df tibble
#'
#'
get_normalized_protein_junc_pos <- function(df){

  df_mod <- df %>%
    # bring all junction positions to left side
    dplyr::mutate(
      normalized_cds_junc_pos = case_when(
        intron_retention & junc_pos_cds_wt == 0 ~
          junc_pos_cds - cds_length_difference,
        cds_length_difference > 0  & junc_pos_cds_wt == 0 ~
          junc_pos_cds - cds_length_difference,
        TRUE ~ junc_pos_cds
      ), normalized_protein_junc_pos = case_when(
        intron_retention & junc_pos_cds_wt == 0 ~
          ceiling(normalized_cds_junc_pos / 3),
        !frame_shift & protein_length_difference > 0  & junc_pos_cds_wt == 0 ~
          protein_junc_pos - protein_length_difference,
        frame_shift & cds_length_difference > 0  & junc_pos_cds_wt == 0 ~
          protein_junc_pos - ceiling(cds_length_difference/3),

        TRUE ~ protein_junc_pos
      ),
      WT_protein_length_difference = ifelse(
        !frame_shift & protein_length_difference > 0 ,
        1,
        abs(protein_length_difference)
      )
    )

  # get left and right WT AA
  df_mod <- df_mod %>%
    dplyr::mutate(
      exon1_end_AA = substr(
        protein,
        normalized_protein_junc_pos,
        normalized_protein_junc_pos
      ),
      exon1_end_AA_WT = substr(
        protein_wt,
        normalized_protein_junc_pos,
        normalized_protein_junc_pos
      ),
      exon2_start_AA_WT = substr(
        protein_wt,
        normalized_protein_junc_pos + WT_protein_length_difference,
        normalized_protein_junc_pos + WT_protein_length_difference
      )
    )

  df_mod <- df_mod %>%
    dplyr::mutate(
      normalized_protein_junc_pos = case_when(
        # left WT-AA
        !is_first_reading_frame & (exon1_end_AA == exon1_end_AA_WT) ~
          ceiling(normalized_cds_junc_pos / 3),
        !is_first_reading_frame & (exon1_end_AA != exon1_end_AA_WT) ~
          floor(normalized_cds_junc_pos / 3),
        intron_retention & (exon1_end_AA != exon1_end_AA_WT) ~
          floor(normalized_cds_junc_pos / 3),
        !is_first_reading_frame  ~
          floor(normalized_cds_junc_pos / 3),
        TRUE ~ normalized_protein_junc_pos
      )
    )

  return(df_mod)
}

#' Get peptide context sequence given a tibble with
#'  `intron_retention`
#'  `normalized_protein_junc_pos`
#'  `is_first_reading_frame`
#'  `frame_shift`
#'  `exon1_end_AA`
#'  `exon1_end_AA_WT`
#'  `exon2_start_AA_WT`,
#'  `protein_len`
#'  `protein_len`
#'  `protein_length_difference`
#'  `protein`
#'
#' @param df Data frame with information about position of junction etc.
#' @param flanking_size number amino acids left and right of the breakpoint or novel sequence part
#'
#'
get_peptide_context <- function(df, flanking_size = 14){

  # peptide start coordinate in full protein sequence
  # 14 AA upstream of junction
  df_mod <- df %>%
    dplyr::mutate(
      pep_start = pmax(normalized_protein_junc_pos - flanking_size + 1, 1)
      ) %>%
    # if not 1st reading frame, there can be an additional novel AA
    mutate(
      addtional_AA = case_when(
       !frame_shift & !is_first_reading_frame &
         (exon1_end_AA != exon1_end_AA_WT & exon1_end_AA != exon2_start_AA_WT )
         ~
         1,
       TRUE ~ 0
      )
    )

  # peptide end coordinate in full protein sequence
  # frame-shifts: full sequence downstream of junction
  # no frame-shift: novel  sequence + WT-AA downstream of junction or until stop codon
  df_mod <- df_mod %>%
    dplyr::mutate(
      pep_end = case_when(
        # IR: WT-AA + ir + WT-AA or until stop
        intron_retention & !frame_shift ~
          pmin(normalized_protein_junc_pos + addtional_AA + flanking_size, protein_len),
        # frameshift: until stop
        frame_shift ~
          protein_len,
        # no frame-shift & no IR & insertion of novel sequence: WT-AA + novel seq + WT-AA or until stop
        !frame_shift & protein_length_difference > 0 ~
          pmin(normalized_protein_junc_pos + addtional_AA + protein_length_difference + flanking_size, protein_len),
        # no frame-shift & no IR & deletion of sequence: WT-AA + junction + WT-AA or until stop
        !frame_shift & protein_length_difference <= 0 ~
          pmin(normalized_protein_junc_pos + addtional_AA + flanking_size, protein_len)
      ))

  # get raw sequence
  df_mod <- df_mod %>%
    dplyr::mutate(peptide_context_seq_raw = stringr::str_sub(protein, start = pep_start, end = pep_end))


  # calculate junction position relative to context sequence
  # get peptide_context sequence until first stop codon
  df_mod <- df_mod %>%
    dplyr::mutate(
      peptide_context_junc_pos =
        normalized_protein_junc_pos - pep_start + 1
    ) %>%
    dplyr::mutate(
      peptide_context =
        seq_truncate_nonstop(peptide_context_seq_raw, peptide_context_junc_pos)
    )

  return(df_mod)
}



