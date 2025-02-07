
#' Annotates splice junctions with resulting CDS and peptide sequence.
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
#' @return A data.frame with the same rows as the input `df` but with the
#'  following additional column(s):
#'
#'  - `cds_mod_id` An identifier made from `tx_id` and `junc_id`
#'  - `junc_pos_cds` the junction position in the modified CDS sequence
#'  - `frame_shift` Indicator whether junction leads to frame shift.
#'  - `is_first_reading_frame` Indicator whether modified CDS sequence is
#'    translated into `protein` sequence using the 1st reading frame
#'    (i.e. reading is starting from the first nucleotide) for in-frame peptides.
#'  - `normalized_cds_junc_pos` The normalized position of the junction in the
#'    modified CDS sequence to the left junction side.
#'  - `protein` The full protein sequence of the translated modified CDS.
#'  - `normalized_protein_junc_pos` The normalized position of the junction in
#'    the `protein` sequence to the left junction side.
#'  - `peptide_context_junc_pos` The junction position relative to the `peptide_context` sequence
#'  - `junc_in_orf` Indicator whether the junction is located in an open reading frame.
#'  - `peptide_context_seq_raw` The peptide sequence around the junction including stop codons.
#'  - `peptide_context` The peptide sequence around the junction truncated after stop codons.
#'  - `truncated_cds` Indicator whether the mutated gene product is a truncated
#'     from of the WT gene product. If TRUE, `peptide_context` = NA.
#'  - `cds_description` Descriptor of of the mutated gene product. Can be one of
#'    c("mutated cds", "truncated wt cds", "no mutated gene product", "no wt cds", "not in ORF")
#'
#'   If the `keep_ranges` is TRUE, the following additional columns are added to
#'   the output data.frame:
#'
#'  - `cds_lst` A list of \code{\link[GenomicRanges]{GRanges}} with
#'        the original CDS as provided in `tx_id` column and `cds` object..
#'  - `cds_mod_lst` A list of \code{\link[GenomicRanges]{GRanges}} with the
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
  is_intron_retention  <- jx@ranges@width == 2

  # modify transcripts by applying the splice junctions
  cds_mod <- modify_tx(cds_lst, jx)

  # get junction position in altered CDS
  junc_pos_cds = get_junc_pos(cds_mod, jx)
  junc_in_cds = !is.na(junc_pos_cds)

  # get junction position in wt CDS; is 0 if
  junc_pos_cds_wt = get_junc_pos(cds_lst, jx)

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
      is_intron_retention = is_intron_retention ,
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
      protein_len = BiocGenerics::width(protein) %>% as.numeric()
    )

  df_positions <- df_positions  %>%
    get_normalized_protein_junc_pos()%>%
    annotate_junc_in_orf()

  # empty mutated proteins to NA --> junc outside of ORF
  # mutated gene product truncated version of WT gene product?
  df_positions <- df_positions %>%
    annotate_truncated_cds()

  df_annotated_peptide <- df_positions %>%
    get_peptide_context(flanking_size = flanking_size)

  # Annotate table
  df_annotated_peptide <- df_annotated_peptide %>%
    dplyr::mutate(
      # add NA for context sequences if the junction position is not in an open reading frame
      # or cds is truncated, i.e. mutated gene product is truncated version of the WT gene product
      peptide_context_seq_raw = ifelse(junc_in_orf & !truncated_cds, as.character(peptide_context_seq_raw), NA),
      peptide_context = ifelse(junc_in_orf & !truncated_cds, as.character(peptide_context), NA),
      peptide_context_junc_pos = ifelse(junc_in_orf & !truncated_cds, peptide_context_junc_pos, NA),
    )

  # if keep_ranges argument is TRUE add list columns of GRanges as transcripts
  if(keep_ranges){
    df_annotated_peptide <- df_annotated_peptide %>%

      dplyr::mutate(
        cds_lst = as.list(cds_lst),
        cds_mod_lst = as.list(cds_mod),
      )
  }

  df_annotated_peptide <- df_annotated_peptide %>%
    dplyr::select(
      -is_intron_retention ,
      -strand,
      -junc_in_cds,
      -protein_len,
      -pep_start,
      -pep_end,
      -exon_end_AA,
      -exon_end_AA_WT,
      -cds_length_difference,
      -additional_novel_AA,
      -junc_pos_cds_wt,
    )

  # add annotations to input data.frame
  df <- df %>%
    dplyr::left_join(df_annotated_peptide, by = c("junc_id", "tx_id")) %>%
    dplyr::mutate(
      cds_description = ifelse(is.na(protein_wt), "no wt cds", cds_description)
    ) %>%
  # return protein seq only until first stop codon
    dplyr::mutate(protein = protein %>% str_extract("[^*]+"))


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
#' @keywords internal
is_first_reading_frame <- function(df){

  df_mod <- df %>%
   dplyr::mutate(
     is_first_reading_frame = case_when(
       round((normalized_cds_junc_pos/3) %% 1, 2) == 0 ~ TRUE,
       round((normalized_cds_junc_pos/3) %% 1, 2) != 0 ~ FALSE,
     )
   )

}


#' Tests if junction position is in ORF
# i.e. no stop codon `*` in whole seq before junction given a tibble with:
#'  `normalized_protein_junc_pos`
#'  `protein_len`
#'  `protein`
#'
#' @param df tibble
#'
#' @keywords internal
annotate_junc_in_orf <- function(df){

  df_mod <- df %>%
    dplyr::mutate(
      protein_until_junction = stringr::str_sub(
        protein,
        start = 1,
        end = pmin(normalized_protein_junc_pos, protein_len)
      ),
      junc_in_orf = stringr::str_detect(protein_until_junction, "\\*", negate = TRUE),
      junc_in_orf = ifelse(is.na(junc_in_orf), FALSE, junc_in_orf)
    ) %>%
    dplyr::select(-protein_until_junction)

  return(df_mod)
}


#' Tests if junction leads to truncated mutated gene product
#   given a tibble with:
#'  `protein_wt`
#'  `protein`
#'  `protein_len`
#'  `junc_in_orf`
#'
#' @param df tibble
#'
#' @keywords internal
annotate_truncated_cds <- function(df){

  df_mod <- df %>%
    # get protein seq until first stop codon
    dplyr::mutate(
      protein = ifelse(protein == "", NA, protein),
      protein_until_stop_codon =  stringr::str_remove(protein, "\\*.*")
    ) %>%
    dplyr::mutate(
      # protein_until_stop_codon cannot be "" or NA as search pattern otherwise warning --> assign "empty"
      protein_until_stop_codon = ifelse(
        protein_until_stop_codon == "" |
          is.na(protein_until_stop_codon),
        "empty",
        protein_until_stop_codon
      ),
      # annotate if protein is truncated, i.e. a wt sequence only
      truncated_cds = stringr::str_detect(fixed(protein_wt), fixed(protein_until_stop_codon)),
      # describe annotated cds
      cds_description = case_when(
        !junc_in_orf ~ "not in ORF",
        protein_until_stop_codon == "empty" ~ "no mutated gene product",
        truncated_cds & cds_length_difference == 0 ~ "wt cds",
        truncated_cds ~ "truncated wt cds",
        TRUE ~ "mutated cds"
      ),
      truncated_cds = ifelse(!junc_in_orf, NA, truncated_cds)
    ) %>%
    select(-protein_until_stop_codon)

  return(df_mod)
}


#' Annotate the normalized junction position in the protein,
# i.e. the postion of the last WT amino acid in the mutated protein on the left side
#'  `frame_shift`
#'  `is_intron_retention`
#'  `protein_junc_pos`
#'  `protein`
#'  `protein_wt`
#'  `junc_pos_cds`
#'  `is_first_reading_frame`
#'  `junc_pos_cds_wt`
#'  `cds_length_difference`
#'
#' @param df tibble
#'
#' @keywords internal
get_normalized_protein_junc_pos <- function(df){

  df_mod <- df %>%
    # bring all junction positions to left side
    # left from normalised junc in cds/protein sequence is only wt
    dplyr::mutate(

      normalized_cds_junc_pos = case_when(
        # IR and junc start is not on cds
        # normalize to junc on left side
        is_intron_retention & junc_pos_cds_wt == 0 ~
          junc_pos_cds - cds_length_difference,

        # Non-IR events (e.g. ASS events) and junc start coordinate is in intron
        # leading to nucleotidesin the cds seq
        # shift junc pos to left to have altered sequence on the right side
        cds_length_difference > 0  & junc_pos_cds_wt == 0 ~
          junc_pos_cds - cds_length_difference,

        # all other cases (e.g. ES and exitrons) should have junc start coordinate
        # on a exon and novel nucleotides in CDS on the left side of the junction
        TRUE ~ junc_pos_cds
      ),

      # First, calculate protein junc pos based on position of junction in cds
      # Depending on th reading frame and the resulting codon this
      # might have a high chance to be a WT AA
      normalized_protein_junc_pos_pre =  ceiling(normalized_cds_junc_pos / 3)
    )%>%
    is_first_reading_frame()


  # get AA at calculated junc position in the altered and wt protein
  df_mod <- df_mod %>%
    dplyr::mutate(
      exon_end_AA = substr(
        protein,
        normalized_protein_junc_pos_pre,
        normalized_protein_junc_pos_pre
      ),
      exon_end_AA_WT = substr(
        protein_wt,
        normalized_protein_junc_pos_pre,
        normalized_protein_junc_pos_pre
      )
    )

  # update normalized_protein_junc_pos based on AA resulting by the codon
  # that is directly affected by the junction if not first reading frame
  df_mod <- df_mod %>%
    dplyr::mutate(

      additional_novel_AA = case_when(
        !is_first_reading_frame & (exon_end_AA != exon_end_AA_WT) ~ 1,
        TRUE ~ 0
      ),

      normalized_protein_junc_pos = case_when(

        # first reading frame
        # use normalized_protein_junc_pos as calculated above
        is_first_reading_frame ~ normalized_protein_junc_pos_pre,

        # no first reading frame but affected codon
        # results in wt aa
        # use normalized_protein_junc_pos as calculated above
        !is_first_reading_frame & additional_novel_AA == 0 ~
          normalized_protein_junc_pos_pre,

        # not first reading frame & affected codon results in a novel AA
        # do not include the novel AA in the position
        !is_first_reading_frame & additional_novel_AA == 1 ~
          floor(normalized_cds_junc_pos / 3),

      )
    )

  return(df_mod)
}

#' Get peptide context sequence given a tibble with
#'  `normalized_protein_junc_pos`
#'  `is_first_reading_frame`
#'  `frame_shift`
#'  `additional_novel_AA`
#'  `protein`
#'  `cds_length_difference`
#'
#' @param df Data frame with information about position of junction etc.
#' @param flanking_size number amino acids left and right of the breakpoint or novel sequence part
#'
#' @keywords internal
get_peptide_context <- function(df, flanking_size = 14){

  # peptide start coordinate in full protein sequence
  # `flanking_size`AA upstream of junction
  df_mod <- df %>%
    dplyr::mutate(
      pep_start = pmax(normalized_protein_junc_pos - flanking_size + 1, 1)
    )

  # peptide end coordinate in full protein sequence
  # frame-shifts: full sequence downstream of junction
  # no frame-shift: novel  sequence + WT-AA downstream of junction or until stop codon
  df_mod <- df_mod %>%
    dplyr::mutate(
      pep_end = case_when(

        # frameshift: until stop
        frame_shift ~ protein_len,

        # no frame-shift & no novel AAs:
        # WT-AA + junction + WT-AA or until stop
        # (e.g. ASS, ES)
        # additional AA if junction changes codon, needs to be considered
        # --> one more WT AA on right sight needed
        !frame_shift & cds_length_difference <= 0 ~
          pmin(normalized_protein_junc_pos + additional_novel_AA + flanking_size, protein_len),

        # no frame-shift  & novel AAs:
        # WT-AA + novel seq + WT-AA or until stop
        # (e.g. ASS, IRs)
        # TODO: do we need additional_novel_AA?
        !frame_shift & cds_length_difference > 0 ~
          pmin(normalized_protein_junc_pos + additional_novel_AA + (cds_length_difference / 3) + flanking_size, protein_len)

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



