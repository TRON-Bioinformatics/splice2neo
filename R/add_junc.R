
#' modifies transcripts by applying splice-junctions
#'
#' @param tx transcripts a named \code{\link[GenomicRanges]{GRangesList}} with
#'    transcripts defined as GRanges of exons.
#' @param jx splice junctions as GRanges objects with the same length as `tx`.
#'
#' @return a \code{\link[GenomicRanges]{GRangesList}} with altered transcripts.
#'
#' The algorithm works in the following steps.
#'
#' - Get transcripts for which the whole transcript range contains `pos1` and
#' `pos2` as GRangesList of GRanges of exons
#' - For each tx in transcripts do
#'   - if `pos1` amd `pos2` are on the same exon, duplicate exon and consider
#'     the first for `pos1` and the second for `pos2`.
#'   - if `pos1` overlaps an exon `e`:
#'      shorten `e` to end at `pos1`
#'    else:
#'      get next exon `e` left of `pos`:
#'      extend `e` to end in `pos1`
#'   - if `pos2` overlaps an exon `e`:
#'      shorten `e` to start at `pos2`
#'    else:
#'      get next exon `e` right of `pos2`:
#'      extend `e` to start in `pos2`
#'   - remove all exons completeley contained in the interval \[`pos1`, `pos2`\].
#'
#' @examples
#'
#' @export
add_junc <- function(tx, jx){

  ## assume same length of tx and jx
  stopifnot(length(tx) == length(jx))

  # get individual GRanges objects for start and end position of junction
  j1 <- GenomicRanges::resize(jx, width = 1, fix="start")
  j2 <- GenomicRanges::resize(jx, width = 1, fix="end")

  # convert GRangesLists as natural list objects
  tx_lst <- as.list(tx)
  jx_lst <- S4Vectors::split(jx, 1:length(jx)) %>% as.list()
  j1_lst <- S4Vectors::split(j1, 1:length(j1)) %>% as.list()
  j2_lst <- S4Vectors::split(j2, 1:length(j2)) %>% as.list()


  # Get indexes of affected exons ==============================================
  # get exon index in transcripts that are affected by the first position of the junction
  exon_idx1 = furrr::future_map2_int(tx_lst, j1_lst, ~GenomicRanges::findOverlaps(.y, .x, select = "first"))

  # get index of next exon left of first junction position
  left_exon_idx = furrr::future_map2_int(tx_lst, j1_lst, ~GenomicRanges::follow(.y, .x, ignore.strand = TRUE))

  # get exon index in transcripts that are affected by the second position of the junction
  exon_idx2 = furrr::future_map2_int(tx_lst, j2_lst, ~GenomicRanges::findOverlaps(.y, .x, select = "first"))

  # get index of next exon right of second junction position
  right_exon_idx = furrr::future_map2_int(tx_lst, j2_lst, ~GenomicRanges::precede(.y, .x, ignore.strand = TRUE))

  # get the indexes of exons that needs to be modified.
  # If position is not on an exon, the index is NA and the left (for junction
  # pos 1) or right (for junction pos 2) exon needs to be considered.
  used_exon_idx1 = ifelse(!is.na(exon_idx1), exon_idx1, left_exon_idx)
  used_exon_idx2 = ifelse(!is.na(exon_idx2), exon_idx2, right_exon_idx)

  # Exitorn: junction within same exon =========================================

  # If both junction coordinates are on the same exon, we need to duplicate the
  # exon in the transcript and modify the start of the first and end of the second
  exitron_junc <- used_exon_idx1 == used_exon_idx2

  # duplicate
  tx_lst[exitron_junc] = furrr::future_map2(tx_lst[exitron_junc], used_exon_idx1[exitron_junc],
                              .f = function(gr, exon_idx){
                                gr <- gr[c(1:exon_idx, exon_idx:length(gr))]
                                return(gr)
                              }
  )

  # for exitron junction, use the second (duplicated) exon to modify the
  # junction end point as new exon start coordinate
  used_exon_idx2[exitron_junc] <- used_exon_idx2[exitron_junc] + 1

  # Modify affected exon start and end coordinates =============================

  # modify exon end position for first junction pos
  tx_lst = furrr::future_pmap(
    list(gr = tx_lst,
         at = used_exon_idx1,
         pos = BiocGenerics::start(jx)),
    .f = function(gr, at, pos){
      BiocGenerics::end(gr[at]) <- pos
      return(gr)
    }
  )

  # modify exon start position for second junction pos
  tx_lst = furrr::future_pmap(
    list(gr = tx_lst,
         at = used_exon_idx2,
         pos = BiocGenerics::end(jx)),
    .f = function(gr, at, pos){
      BiocGenerics::start(gr[at]) <- pos
      return(gr)
    }
  )

  # remove all exons in between ================================================

  # Build rangs of *intronic* space between junctions, +/-1bp inwards of the exon start/end postions
  jx_space_lst <- jx %>%
    IRanges::narrow(start = 2, end = -2) %>%

    # convert into list
    S4Vectors::split(1:length(.)) %>% as.list()

  # iterate over transcripts and keep only exons that do not be contained in the junction range
  tx_lst <- furrr::future_map2(tx_lst, jx_space_lst, ~.x[!IRanges::overlapsAny(.x, .y, type="within")])

  tx_out <- GenomicRanges::GRangesList(tx_lst)

  return(tx_out)

}

#' update end coordinate of a subset of ranges by given indices
#'
#' @param grl A \code{\link[GenomicRanges]{GRangesList}} object
#' @param pos The new coordinate
#' @param at The indices of ranges to update
#'
#' @return A modified list object
grl_update_end_at <- function(grl, at, pos){

  stopifnot(length(grl) == length(at) | length(at) == 1)
  stopifnot(length(grl) == length(pos) | length(pos) == 1)

  furrr::future_pmap(
    list(gr = as.list(grl),
         at = at,
         pos = pos),
    .f = function(gr, at, pos){
            BiocGenerics::end(gr[at]) <- pos
            return(gr)
          }
  )
}


#' update start coordinate of a subset of ranges by given indices
#'
#' @param grl A \code{\link[GenomicRanges]{GRangesList}} object
#' @param pos The new coordinate
#' @param at The indices of ranges to update
#'
#' @return A modified list object
grl_update_start_at <- function(grl, at, pos){

  stopifnot(length(grl) == length(at) | length(at) == 1)
  stopifnot(length(grl) == length(pos) | length(pos) == 1)

  furrr::future_pmap(
    list(gr = as.list(grl),
         at = at,
         pos = pos),
    .f = function(gr, at, pos){
      BiocGenerics::start(gr[at]) <- pos
      return(gr)
    }
  )
}

#' adds junction position in transcripts
#'
#' @param tx transcripts a \code{\link[GenomicRanges]{GRangesList}} with
#' transcripts defined as GRanges of exons.
#' @param jx splice junctions as GRanges objects
#'
#' @return an integer vector with junction positions
#'
#' @export
add_junc_pos <- function(tx, jx){
  # TODO
}
