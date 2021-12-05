
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

  # build common seqlevels
  common_seqlevels <- unique(c(
    GenomeInfoDb::seqlevels(j_start),
    GenomeInfoDb::seqlevels(tx),
    "Z" # add custom seqlevel for zero-size GRangs later
  ))

  GenomeInfoDb::seqlevels(j_start) <- common_seqlevels
  GenomeInfoDb::seqlevels(tx) <- common_seqlevels

  # check if junction overlap the whole range of the transcript ---------------

  # buid full range of transcript
  range_tx <- base::range(tx)

  # check which transcripts have zero length, i.e. consists only of empty ranges
  # these empty ranges get lost when converting to GRanges obejects and therefore
  # will be replaced by an zeor-size toy range object
  empty_ranges <- range_tx %>% GenomicRanges::width() %>% as.integer() %>% is.na()
  if (any(empty_ranges)){
    range_tx[which(empty_ranges)] <-  GenomicRanges::GRangesList(
      GenomicRanges::GRanges("Z", IRanges::IRanges(0, 0), "+")
    )
  }

  range_tx <-  unlist(range_tx)

  on_tx <- IRanges::poverlaps(j_start, range_tx) %>%
    as.logical()

  # initialize with NA
  pos_tx <- rep(NA, length(tx))

  # map junction positions on transcript sequence position
  suppressWarnings(
    pos_tx[on_tx] <- GenomicFeatures::pmapToTranscripts(j_start[on_tx], tx[on_tx]) %>%
      BiocGenerics::start()
  )

  return(pos_tx)
}
