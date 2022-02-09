
#' Parse .csv file from MMsplice output as data.frame
#'
#' @param infile path to a .csv file, the output of MMsplice
#' @return a [tibble][tibble::tibble-package] with one row per variant and
#' affected exon.
#'
#' @examples
#'
#' in_file <- system.file("extdata", "mmsplice_pred.csv", package = "splice2neo")
#' parse_mmsplice(in_file)
#'
#' @export
parse_mmsplice <- function(infile){

  # infile <- "https://raw.githubusercontent.com/gagneurlab/MMSplice_MTSplice/master/tests/data/pred.csv"

  df <- readr::read_csv(infile, col_types = readr::cols(
    .default = readr::col_double(),
    ID = readr::col_character(),
    exons = readr::col_character(),
    exon_id = readr::col_character(),
    gene_id = readr::col_character(),
    gene_name = readr::col_character(),
    transcript_id = readr::col_character()
  ))

  # add mut_id
  df <- df %>%
    mutate(
      mut_id = str_split(ID, ":|>", 4) %>%
        furrr::future_map_chr(~str_c(.x[1], "_", .x[2], "_", .x[3], "_", .x[4]))
    )

  return(df)
}


#' Annotates the mmsplice output with additional columns including the junction as junc_id.
#'
#' @param mmsplice_df A data.frame like object from mmsplice output. This should at least
#'   have the following columns:
#'     - ID
#'     - exon_id
#'     - delta_logit_psi
#'
#' @param transcripts a GRangesList with transcripts defined as GRanges of exons
#'   created by `GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)`.
#'
#' @return A data.frame like object like the input but with additional columns:
#' - junc_id: `<chr>_<pos1>_<pos2>_<strand>`
#' - effect: either `exon_skipping` or `exon_inclusion`
#'
#'
#' If the logit delta PSI score is <= 0, the event is treated as exon skipping.
#' In this case a junction is build from the end of the exon before, and the
#' start of the exon after that.
#'
#' If the logit delta PSI score is > 0, the event is treated as cassette exon.
#' In this case a (canonical) junction is build from the end of the upstream
#' exon to the start of the affected exon and from the end of the affected exon
#' to the start of the next exon.
#'
#'
#' MMsplice predicts the change on percent spliced in (PSI) for a given
#' annotated exon. Therefore, only exon inclusion and exon skipping needs to be
#' converted to junctions.
#'
#'  - For *exon skipping* the end of the previous and start of the next exon
#' build the junction
#'  - For *exon inclusion* the (canonical) junction is build
#' from the end of the upstream exon to the start of the affected exon and from
#' the end of the affected exon  to the start of the next exon.
#'
#'
#' @export
annotate_mmsplice <- function(mmsplice_df, transcripts){

  mmsplice_annot <- mmsplice_df %>%

    # select(nr, mut_id, ID, exons, exon_id, gene_id, transcript_id, delta_logit_psi) %>%

    # filter out a few transcripts that are not in the txdb object
    filter(transcript_id %in% names(transcripts)) %>%

    mutate(
      effect = ifelse(delta_logit_psi <= 0, "exon_skipping", "exon_inclusion"),

      # get junctions for exon skipping
      skip_junc = get_exon_skipping_junction(exon_id, transcript_id, transcripts),

      # get junctions for exon inclusion
      incl_junc = get_exon_inclusion_junction(exon_id, transcript_id, transcripts),

      # depending on effect take proper junctions
      junc_id_lst = case_when(
        effect == "exon_skipping" ~ skip_junc,
        effect == "exon_inclusion" ~ incl_junc,
      )
    )

  # tictoc::toc()
  #

  mmsplice_junc_df <- mmsplice_annot %>%
    select(-skip_junc, -incl_junc) %>%
    unnest(junc_id_lst) %>%
    rename(junc_id = junc_id_lst)%>%
    filter(!is.na(junc_id))

  return(mmsplice_junc_df)

}


#' Compute the resulting junction (junc_id) for exon skipping of given exon
#' and transcript.
#'
#' @param exon_id A vector of exon IDs
#' @param transcript_id A vector of transcript IDs
#' @param transcripts a GRangesList with transcripts defined as GRanges of exons
#'   created by `GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)`.
#'
#' @return A list of junction IDs formatted as `<chr>_<left>_<right>_<strand>`.
#'
#' @examples
#'
#'transcript_id <- c("ENST00000243347", "ENST00000460812")
#'exon_id <- c("ENSE00000840477", "ENSE00003481758")
#'transcripts <- toy_transcripts
#'
#'get_exon_skipping_junction(exon_id, transcript_id, transcripts)
#'
#' @export
get_exon_skipping_junction <- function(exon_id, transcript_id, transcripts){

  # get subset of transcripts
  tx_sub <- transcripts[transcript_id]

  # get lsit of exon IDs
  # toy_transcripts %>% BiocGenerics::lapply(function(x){x$exon_name = x$exon_name %>% str_extract("^ENSE\\d{7}")})
  exon_lst <- tx_sub %>% BiocGenerics::lapply(function(x){x$exon_name})

  # get index of exon in transcript
  exon_idx <- purrr::map2_int(exon_id, exon_lst, match)

  # get total exon number per transcript
  exon_n <- purrr::map_int(exon_lst, length)

  # get strand of transcript to select proper left and right exon
  exon_strand = BiocGenerics::strand(tx_sub) %>% sapply(unique) %>% as.character()
  strand_offset <- ifelse(exon_strand == "+", 1, -1)

  # get left and right exon index
  exid_left_idx <- exon_idx - strand_offset
  exid_left_idx <- ifelse(exid_left_idx >= 1 & exid_left_idx <= exon_n, exid_left_idx, NA)
  exid_right_idx <- exon_idx + strand_offset
  exid_right_idx <- ifelse(exid_right_idx >= 1 & exid_right_idx <= exon_n, exid_right_idx, NA)

  # get start and end coordinates of all exons
  exon_ends <- S4Vectors::end(tx_sub) %>% as.list()
  exon_starts <- S4Vectors::start(tx_sub) %>% as.list()

  # get end of left exon and start of right exon
  left_end <- purrr::map2_int(exon_ends, exid_left_idx, ~.x[.y])
  right_start <- purrr::map2_int(exon_starts, exid_right_idx, ~.x[.y])

  # get strand and chromosome for all exons
  # exon_strand = purrr::map2_chr(as.list(strand(tx_sub)), exon_idx, ~as.character(.x[.y]))
  exon_chr = purrr::map2_chr(as.list(GenomeInfoDb::seqnames(tx_sub)), exon_idx, ~ifelse(is.na(.y), NA, as.character(.x[.y])))

  junc_id <- generate_junction_id(exon_chr, left_end, right_start, exon_strand)

  # remove NA's
  # junc_id <- junc_id[!is.na(junc_id)]

  return(as.list(junc_id))
}

#' Compute the resulting junctions (junc_id) for exon inclusion of given exon
#' and transcript.
#'
#' @param exon_id A vector of exon IDs
#' @param transcript_id A vector of transcript IDs
#' @param transcripts a GRangesList with transcripts defined as GRanges of exons
#'   created by `GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)`.
#'
#' @return A list of junction IDs formatted as `<chr>_<left>_<right>_<strand>`.
#'
#' @examples
#'
#'transcript_id <- c("ENST00000243347", "ENST00000460812")
#'exon_id <- c("ENSE00000840477", "ENSE00003481758")
#'transcripts <- toy_transcripts
#'
#'get_exon_inclusion_junction(exon_id, transcript_id, transcripts)
#'
#' @export
get_exon_inclusion_junction <- function(exon_id, transcript_id, transcripts){

  # get subset of transcripts
  tx_sub <- transcripts[transcript_id]

  # get lsit of exon IDs
  exon_lst <- tx_sub %>% BiocGenerics::lapply(function(x){x$exon_name})

  # get index of exon in transcript
  exon_idx <- purrr::map2_int(exon_id, exon_lst, match)

  # get total exon number per transcript
  exon_n <- purrr::map_int(exon_lst, length)

  # get strand of transcript to select proper left and right exon
  exon_strand = BiocGenerics::strand(tx_sub) %>% sapply(unique) %>% as.character()
  strand_offset <- ifelse(exon_strand == "+", 1, -1)

  # get left and right exon index
  exid_left_idx <- exon_idx - strand_offset
  exid_left_idx <- ifelse(exid_left_idx >= 1 & exid_left_idx <= exon_n, exid_left_idx, NA)
  exid_right_idx <- exon_idx + strand_offset
  exid_right_idx <- ifelse(exid_right_idx >= 1 & exid_right_idx <= exon_n, exid_right_idx, NA)

  # get start and end coordinates of all exons
  exon_ends <- S4Vectors::end(tx_sub) %>% as.list()
  exon_starts <- S4Vectors::start(tx_sub) %>% as.list()

  # get end of left exon and start of right exon
  left_end <- purrr::map2_int(exon_ends, exid_left_idx, ~.x[.y])
  exon_start <- purrr::map2_int(exon_starts, exon_idx, ~.x[.y])

  exon_end <- purrr::map2_int(exon_ends, exon_idx, ~.x[.y])
  right_start <- purrr::map2_int(exon_starts, exid_right_idx, ~.x[.y])

  # get strand and chromosome for all exons
  # exon_strand = purrr::map2_chr(as.list(strand(tx_sub)), exon_idx, ~as.character(.x[.y]))
  # exon_chr = purrr::map2_chr(as.list(seqnames(tx_sub)), exon_idx, ~as.character(.x[.y]))
  # exon_strand = purrr::map2_chr(as.list(strand(tx_sub)), exon_idx, ~as.character(.x[.y]))
  exon_chr = purrr::map2_chr(as.list(GenomeInfoDb::seqnames(tx_sub)), exon_idx, ~ifelse(is.na(.y), NA, as.character(.x[.y])))

  junc_left <- generate_junction_id(exon_chr, left_end, exon_start, exon_strand)
  junc_right <- generate_junction_id(exon_chr, exon_end, right_start, exon_strand)

  junc_id_lst <- purrr::map2(junc_left, junc_right, c)
  return(junc_id_lst)
}

