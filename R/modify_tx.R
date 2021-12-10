

#' Modify transcript by introducing splice junctions
#'
#' @param tx transcripts a named \code{\link[GenomicRanges]{GRangesList}} with
#'    transcripts defined as GRanges of exons.
#' @param jx splice junctions as GRanges objects with the same length as `tx`.
#'
#' @return a \code{\link[GenomicRanges]{GRangesList}} with altered transcripts.
#'
#' @examples
#' tx <- GenomicRanges::GRangesList(
#'   list(GenomicRanges::GRanges(c(
#'      "1:2-3:+",
#'      "1:5-6:+",
#'      "1:10-15:+"
#'   )))
#' )
#'
#' jx <- GenomicRanges::GRanges(c("1:3-10:+"))
#'
#' modify_tx(tx, jx)
#'
#' @export
modify_tx <- function(tx, jx){

  ## assume same length of tx and jx
  stopifnot(length(tx) == length(jx))

  # convert to GRangesList
  tx <- GenomicRanges::GRangesList(tx)

  # bring transcripts and junction to common seqlevels
  seq_levels <- union(GenomeInfoDb::seqlevels(tx), GenomeInfoDb::seqlevels(jx))
  GenomeInfoDb::seqlevels(jx) <- seq_levels
  GenomeInfoDb::seqlevels(tx) <- seq_levels

  # convert GRangesLists as natural list objects
  tx_lst <- as.list(tx)
  jx_grl <- S4Vectors::split(jx, 1:length(jx))
  jx_lst <- jx_grl %>% as.list()

  # get full range of each transcript as GRanges object
  tx_range <- unlist(range(tx), use.names=FALSE)

  # get the introns of transcripts as GRangesList
  # this is suggested here:  https://github.com/Bioconductor/GenomicRanges/issues/17
  int <- GenomicRanges::psetdiff(tx_range, tx)
  int_lst <- as.list(int)

  # get introns that overlap with jx
  old_int <- furrr::future_map2(int_lst, jx_lst, IRanges::subsetByOverlaps) %>%
    GenomicRanges::GRangesList()

  # remove overlapping introns
  int <- GenomicRanges::setdiff(int, old_int)

  # resize junction to 1-based explicit intron coordinates
  jx_grl <- jx_grl %>%
    GenomicRanges::resize(IRanges::width(jx_grl) - 2, fix = "center")

  # add jx as intron
  int <- GenomicRanges::union(int, jx_grl)

  # convert back to exons
  exons <- GenomicRanges::psetdiff(tx_range, int)

  # reorder exon ranks for transcripts on negative strand
  exons <- S4Vectors::revElements(exons, any(BiocGenerics::strand(exons) == "-"))

  return(exons)
}
