
#' Annotate splice junctions with affected transcript IDs
#'
#' @param df A data.frame with splice junctions in rows and at least the colums:
#'
#'   -  `junc_id` junction id consisting of genomic coordinates
#'
#' @param transcripts \code{\link[GenomicRanges]{GRangesList}} of transcripts
#'
#' @return A data.frame as the input but with potentially multiple rows
#'     and with the additional column(s):
#'
#'  - `tx_id` the ID of the affected transcript
#'  - `tx_lst` a list of \code{\link[GenomicRanges]{GRanges}} with the transcript
#'
#' @examples
#'
#' @import dplyr
#' @export
add_tx <- function(df, transcripts){

  stopifnot(is.data.frame(df))
  stopifnot("junc_id" %in% names(df))

  # build GRanges of junction
  jx <- junc_to_gr(df$junc_id)

  # full ranges of all transcript as a single GRanges object
  tx_ranges <- base::range(transcripts)

  # Find transcript ranges that contain the junction range
  hits <- GenomicRanges::findOverlaps(jx, tx_ranges, type = "within")
  jx_idx <- hits %>% S4Vectors::queryHits()
  tx_idx <- hits %>% S4Vectors::subjectHits()

  # build data.frame to associate junc_id to transcripts
  junc_to_tx <- tibble::tibble(
    junc_id = df$junc_id[jx_idx],
    tx_id = names(transcripts)[tx_idx],
    tx_lst = as.list(transcripts[tx_idx])
  )

  # join the original input data.frame with the association of junction to transcripts
  out_df <- df %>%
    dplyr::left_join(junc_to_tx, by = "junc_id")

  return(out_df)
}