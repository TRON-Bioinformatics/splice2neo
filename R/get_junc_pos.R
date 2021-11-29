
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

  # check if junction overlap the whole range of the transcript
  range_tx <- unlist(base::range(tx))
  seqlevels(j_start) <- seqlevels(range_tx)
  on_tx <- IRanges::poverlaps(j_start,range_tx ) %>%
    as.logical()

  # initialize with NA
  pos_tx <- rep(NA, length(tx))

  # map junction positions on transcript sequence position
  suppressWarnings(
    pos_tx[on_tx] <- GenomicFeatures::pmapToTranscripts(j_start[on_tx], tx[on_tx]) %>%
      BiocGenerics::start()
  )
  # # if j_start does not map to the transcript 0 is returned (see ?GenomicFeatures::pmapToTranscripts)
  # # This needs to be replaced by NA
  # pos_tx <- pos_tx %>%
  #   dplyr::na_if(0)

  return(pos_tx)
}
