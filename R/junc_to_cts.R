
#' Annotate splice junctions with resulting transcripts and peptide sequence
#'
#' @param junc_id junction id consisting of genomic coordinates
#' @param transcripts \code{\link[GenomicRanges]{GRangesList}} of transcripts
#' @param tx_id a vector of transcript ids in the same length as `junc_id`. If
#'    provided, only transcripts with the associated id will be considered for
#'    each junction. If NA (default) no input transcripts will be considered.
#' @param size the total size of the output sequence (might be shorter if transcripts is shorter)
#' @param bsg \code{\link[BSgenome]{BSgenome}} object such as
#'  \code{\link[BSgenome.Hsapiens.UCSC.hg19]{BSgenome.Hsapiens.UCSC.hg19}}
#'
#' @return A data.frame with annotated input junctions with colums:
#'
#'  - `junc_id` the input `junc_id`
#'  - `tx_id` the id of the used transcript
#'  - `junc_pos_tx` the junction position in the modified transcript sequence
#'  - `tx_id_alt` an identifier made from `tx_id` and `junc_id`
#'  - `cts_seq` the context sequence
#'  - `cts_junc_pos` the junction position in the context sequence
#'  - `cts_junc_pos` the size of the context sequence
#'  - `cts_id` a unique id for the context sequence as hash value using the
#'   XXH128 hash algorithm
#'
#' @examples
#'
#' requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
#' bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
#'
#' junc_to_cts(toy_junc_id, toy_transcripts, tx_id = toy_junc_id_enst,
#'    size = 400, bsg = bsg)
#'
#'@import dplyr
#'
#' @export
junc_to_cts <- function(junc_id, transcripts, tx_id = NA, size = 400, bsg = NULL){

  if(is.null(bsg)){
    message("INFO: Use default genome sequence from BSgenome.Hsapiens.UCSC.hg19")
    bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  }

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

  # compute affected transcripts
  transcript_df <- junc_df %>%
    # filter out NA junctions
    filter(!is.na(junc_id)) %>%
    mutate(
      gr = furrr::future_map_if(.x = tx_id_input,
                                .p = ~!is.na(.x),
                                .f = ~transcripts[.x],
                                .else = function(){transcripts}
                                ),
      tx_df = furrr::future_pmap(
        list(chr, pos1, pos2, gr),
        junc_to_tx)
    )

  # build junction and transcripts specific data sets
  cont_df <- transcript_df %>%
    unnest(tx_df) %>%
    mutate(
      tx_id_alt = str_c(tx_id, "|", junc_id)
    )

  if(nrow(cont_df) > 0){
    peptide_junc_pos <- ceiling(cont_df$junc_pos_tx / 3)
  }else{
    peptide_junc_pos <- numeric()
  }

  # prepare transcripts and junctoin position depending on 0 or more input data
  if(nrow(cont_df) > 0){
    tx_grl <- GenomicRanges::GRangesList(cont_df$tx)
    junc_pos_tx <- cont_df$junc_pos_tx
  }else{
    tx_grl <- GenomicRanges::GRangesList()
    junc_pos_tx <- numeric()
  }

  # get transcript sequence and junction position in sequence
  tx_seq <- GenomicFeatures::extractTranscriptSeqs(bsg, tx_grl)

  # Get context-sequence around junction
  tx_len <- BiocGenerics::width(tx_seq)
  cts_start <- pmax(junc_pos_tx - (size/2) + 1, 1)
  cts_end <- pmin(junc_pos_tx + (size/2), tx_len)

  # extract context sequence from full transcript
  cts_seq <- XVector::subseq(tx_seq, start = cts_start, end = cts_end)

  # calculate junction position relative to context sequence
  cts_junc_pos <- junc_pos_tx - cts_start

  # Annotate table
  cont_df <- cont_df %>%
    mutate(
      cts_seq = as.character(cts_seq),
      cts_junc_pos = cts_junc_pos,
      cts_size = BiocGenerics::width(cts_seq),
      cts_id = furrr::future_map_chr(cts_seq, rlang::hash)
    )

  # join with original input data.frame to preserve rows with NA in `junc_id`
  junc_df %>%
    left_join(
      cont_df,
      by = c("junc_id", "chr", "pos1", "pos2", "strand", "tx_id_input")
    ) %>%
    select(junc_id, tx_id, junc_pos_tx, tx_id_alt, cts_seq, cts_junc_pos,
           cts_size, cts_id)

}



