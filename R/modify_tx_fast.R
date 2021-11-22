

#' Modify transcript by introducing splice junctions
#'
#' @param tx transcripts a named \code{\link[GenomicRanges]{GRangesList}} with
#'    transcripts defined as GRanges of exons.
#' @param jx splice junctions as GRanges objects with the same length as `tx`.
#'
#' @return a \code{\link[GenomicRanges]{GRangesList}} with altered transcripts.
#'
#'
#' TODO:
#' - [x] resize jx to intron (+/-1)
#' - [ ] test for intron retention
#'
modify_tx_fast <- function(tx, jx){

  # convert GRangesLists as natural list objects
  tx_lst <- as.list(tx)
  jx_grl <- S4Vectors::split(jx, 1:length(jx))
  jx_lst <- jx_grl %>% as.list()

  # https://github.com/Bioconductor/GenomicRanges/issues/17

  # get full range of each transcript as GRanges object
  tx_range <- unlist(range(tx), use.names=FALSE)

  # get the introns of transcripts as GRangesList
  int <- GenomicRanges::psetdiff(tx_range, tx)
  int_lst <- as.list(int)


  # get introns that overlap with jx
  old_int <- purrr::map2(int_lst, jx_lst, IRanges::subsetByOverlaps) %>%
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

}

# mod_tx(tx, jx)
# introns <- psetdiff(unlist(range(grl),use.names=FALSE),grl)
# #
# mergen overlapping introns for all tx
