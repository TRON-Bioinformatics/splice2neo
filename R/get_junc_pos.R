
#' Get the position of input junction in the transcript sequences
#'
#' @param tx transcripts a \code{\link[GenomicRanges]{GRangesList}} with
#' transcripts defined as GRanges of exons.
#' @param jx splice junctions as GRanges objects
#'
#' @return an integer vector with junction positions
#'
#' @export
get_junc_pos <- function(tx, jx){

  # assume same length of tx and jx
  stopifnot(length(tx) == length(jx))

  # convert to GRangesList
  if(! class(tx) %in% c("GRangesList", "CompressedGRangesList")){
    tx <- GenomicRanges::GRangesList(tx)
  }

  # get individual GRanges objects for start and end position of junction
  j_start <- GenomicRanges::resize(jx, width = 1, fix="start")

  # map junction positions on transcript sequence position
  pos_tx <- GenomicFeatures::pmapToTranscripts(j_start, tx) %>%
    BiocGenerics::start()

  return(pos_tx)
}
