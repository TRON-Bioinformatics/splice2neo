
#' Annotate splice junctions with resulting transcript sequence
#'
#' @param df A data.frame with splice junctions in rows and at least the columns:
#'
#'   -  `junc_id` junction id consisting of genomic coordinates
#'   -  `tx_id` the ID of the affected transcript (see \code{\link{add_tx}})
#'
#' @param transcripts  as a named \code{\link[GenomicRanges]{GRangesList}} of transcripts
#' @param size the size of the output sequence around the junction position (might be shorter if transcripts is shorter)
#' @param bsg \code{\link[BSgenome]{BSgenome}} object such as
#'  \code{\link[BSgenome.Hsapiens.UCSC.hg19]{BSgenome.Hsapiens.UCSC.hg19}}
#' @param keep_ranges Should GRanges of transcripts and modified transcript be
#' kept? If TRUE, the list columns `tx_lst` and `tx_mod_lst` are added to the output.
#'
#' @return A data.frame with the same rows as the input `df` but with the
#'  following additional column(s):
#'
#'  - `tx_mod_id` an identifier made from `tx_id` and `junc_id`
#'  - `junc_pos_tx` the junction position in the modified transcript sequence
#'  - `cts_seq` the context sequence
#'  - `cts_junc_pos` the junction position in the context sequence
#'  - `cts_size` the size of the context sequence
#'  - `cts_id` a unique id for the context sequence as hash value using the
#'   XXH128 hash algorithm
#'
#'   If the `keep_ranges` is TRUE, the following additional columns are added to
#'   the output data.frame:
#'
#'  - `tx_lst` a list of \code{\link[GenomicRanges]{GRanges}} with
#'        the original transcript as provided in `tx_id` column and `transcripts` object..
#'  - `tx_mod_lst` a list of \code{\link[GenomicRanges]{GRanges}} with
#'        the modified transcript (see \code{\link{modify_tx}})
#'
#' @examples
#'
#' requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
#' bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#'
#' add_context_seq(toy_junc_df, toy_transcripts, size = 20, bsg = bsg)
#'
#' @import dplyr
#'
#' @export
add_context_seq <- function(df,
                            transcripts,
                            size = 400,
                            bsg = NULL,
                            keep_ranges = FALSE){

  stopifnot(is.data.frame(df))
  stopifnot("junc_id" %in% names(df))
  stopifnot("tx_id" %in% names(df))
  stopifnot(class(transcripts) %in% c("GRangesList", "CompressedGRangesList"))
  stopifnot(is.logical(keep_ranges) & length(keep_ranges) == 1)

  # take only the columns junc_id and tx_id and build unique combinations
  df_sub <- df %>%
    dplyr::distinct(junc_id, tx_id) %>%
    # filter for subset with matching tx_id in input transcripts
    dplyr::filter(tx_id %in% names(transcripts))

  # get GRanges as subset of transcripts
  tx_lst <- transcripts[df_sub$tx_id]

  if(is.null(bsg)){
    message("INFO: Use default genome sequence from BSgenome.Hsapiens.UCSC.hg19")
    bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  }

  # get junctions as GRanges object
  jx <- junc_to_gr(df_sub$junc_id)

  # test if junction is an intron retention event
  # TODO: are non IR junctions possible that follow the rule chr:pos-pos+1:strand?
  intron_retention <- ifelse(jx@ranges@width == 2, TRUE, FALSE)

  # modify transcripts by appliing the splice junctions
  tx_mod <- modify_tx(tx_lst, jx)

  # get junction position in altered transcript
  junc_pos_tx = get_junc_pos(tx_mod, jx)

  # only relevant for intron retention events: get the other position of retained
  # intervall
  junc_end_tx <- get_intronretention_alt_pos(tx_lst, tx_mod, jx, intron_retention)

  # get transcript sequence and junction position in sequence
  tx_seq <- GenomicFeatures::extractTranscriptSeqs(bsg, tx_mod)

  # Get context-sequence around junction
  tx_len <- BiocGenerics::width(tx_seq)
  cts_start <- ifelse(intron_retention & junc_end_tx < junc_pos_tx,
                      pmax(junc_end_tx - (size/2) + 1, 1),
                      pmax(junc_pos_tx - (size/2) + 1, 1))
  cts_end <- ifelse(intron_retention & junc_end_tx > junc_pos_tx,
                    pmin((junc_end_tx - 1) + (size/2), tx_len) ,
                    pmin(junc_pos_tx + (size/2), tx_len))


  # if subset positions are NA, set the enrie sequence to NA to force NA
  # in the ouptut of the folllwoing call to XVector::subseq
  tx_seq[is.na(cts_start) | is.na(cts_end)] <- as("", "DNAStringSet")

  # extract context sequence from full transcript
  cts_seq <- XVector::subseq(tx_seq, start = cts_start, end = cts_end) %>%
    as.character() %>%
    dplyr::na_if("") %>%
    as.character()

  # calculate junction position relative to context sequence
  cts_junc_pos <- ifelse(intron_retention & junc_end_tx < junc_pos_tx,
                         junc_end_tx - cts_start + 1,
                         junc_pos_tx - cts_start + 1)
  #cts_junc_pos <- junc_pos_tx - cts_start + 1
  cts_end_pos <- ifelse(intron_retention & junc_end_tx < junc_pos_tx,
         junc_pos_tx - cts_start + 1,
         junc_end_tx - cts_start)

  # only relevant for IRs: get intervals for intron retentions: start_cts,start_IR,end_IR,end_cts
  cts_intervall <- do.call(paste, list(0, cts_junc_pos, cts_end_pos, nchar(cts_seq), sep = ","))

  # Annotate table
  df_sub <- df_sub %>%
    dplyr::mutate(
      tx_mod_id = stringr::str_c(tx_id, "|", junc_id),
      junc_pos_tx = junc_pos_tx,
      cts_seq = as.character(cts_seq),
      cts_junc_pos = ifelse(intron_retention, cts_intervall, cts_junc_pos),
      cts_size = stringr::str_length(cts_seq),
      cts_id = purrr::map2_chr(cts_seq, cts_junc_pos, ~rlang::hash(c(.x, .y)))
    )

  # if keep_ranges argument is TRUE add list columns of GRanges as transcripts
  if(keep_ranges){
    df_sub <- df_sub %>%
      dplyr::mutate(
        tx_lst = as.list(tx_lst),
        tx_mod_lst = as.list(tx_mod),
      )
  }

  # add annotations to input data.frame
  df <- df %>%
    dplyr::left_join(df_sub, by = c("junc_id", "tx_id"))

  return(df)

}
