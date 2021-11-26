
#' Annotate splice junctions with resulting transcript sequence
#'
#' @param df A data.frame with splice junctions in rows and at least the columns:
#'
#'   -  `junc_id` junction id consisting of genomic coordinates
#'   -  `tx_id` the ID of the affected transcript (see \code{\link{add_tx}})
#'   -  `tx_lst` a list of \code{\link[GenomicRanges]{GRanges}} with the
#'         transcript (see \code{\link{add_tx}})
#'
#' @param size the size of the output sequence around the junction position (might be shorter if transcripts is shorter)
#' @param bsg \code{\link[BSgenome]{BSgenome}} object such as
#'  \code{\link[BSgenome.Hsapiens.UCSC.hg19]{BSgenome.Hsapiens.UCSC.hg19}}
#'
#' @return A data.frame as the input with the additional column(s):
#'
#'   - `tx_alt_lst` a list of \code{\link[GenomicRanges]{GRanges}} with
#'        the modified transcript (see \code{\link{modify_tx}})
#'  - `tx_id_alt` an identifier made from `tx_id` and `junc_id`
#'  - `junc_pos_tx` the junction position in the modified transcript sequence
#'  - `cts_seq` the context sequence
#'  - `cts_junc_pos` the junction position in the context sequence
#'  - `cts_size` the size of the context sequence
#'  - `cts_id` a unique id for the context sequence as hash value using the
#'   XXH128 hash algorithm
#'
#' @examples
#'
#' requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
#' bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#'
#' junc_df <- toy_junc_df %>%
#'   dplyr::mutate(
#'     tx_lst = as.list(toy_transcripts[tx_id])
#'   )
#'
#' add_context_seq(junc_df, size = 30, bsg = bsg)
#'
#' @import dplyr
#'
#' @export
add_context_seq <- function(df, size = 400, bsg = NULL){

  stopifnot(is.data.frame(df))
  stopifnot("junc_id" %in% names(df))
  stopifnot("tx_id" %in% names(df))
  stopifnot("tx_lst" %in% names(df))

  if(is.null(bsg)){
    message("INFO: Use default genome sequence from BSgenome.Hsapiens.UCSC.hg19")
    bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  }

  # get junctions as GRanges object
  jx <- junc_to_gr(df$junc_id)

  # modify transcripts by appliing the splice junctions
  tx_alt <- modify_tx(df$tx_lst, jx)

  # build new id
  tx_id_alt = str_c(df$tx_id, "|", df$junc_id)

  # get junction position in altered transcript
  junc_pos_tx = get_junc_pos(tx_alt, jx)

  # get transcript sequence and junction position in sequence
  tx_seq <- GenomicFeatures::extractTranscriptSeqs(bsg, tx_alt)

  # Get context-sequence around junction
  tx_len <- BiocGenerics::width(tx_seq)
  cts_start <- pmax(junc_pos_tx - (size/2) + 1, 1)
  cts_end <- pmin(junc_pos_tx + (size/2), tx_len)

  # if subset positions are NA, set the enrie sequence to NA to force NA
  # in the ouptut of the folllwoing call to XVector::subseq
  tx_seq[is.na(cts_start) | is.na(cts_end)] <- NA_character_

  # extract context sequence from full transcript
  cts_seq <- XVector::subseq(tx_seq, start = cts_start, end = cts_end) %>%
    as.character() %>%
    dplyr::na_if("") %>%
    as.character()

  # calculate junction position relative to context sequence
  cts_junc_pos <- junc_pos_tx - cts_start

  # Annotate table
  df %>%
    dplyr::mutate(
      tx_alt_lst = as.list(tx_alt),
      tx_id_alt = tx_id_alt,
      junc_pos_tx = junc_pos_tx,
      cts_seq = as.character(cts_seq),
      cts_junc_pos = cts_junc_pos,
      cts_size = stringr::str_length(cts_seq),
      cts_id = furrr::future_map_chr(cts_seq, rlang::hash)
    )

}



