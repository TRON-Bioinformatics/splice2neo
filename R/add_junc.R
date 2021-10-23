
#' modifies transcripts by applying splice-junctions
#'
#'
#' @param tx transcripts a named \code{\link[GenomicRanges]{GRangesList}} with
#' transcripts defined as GRanges of exons.
#' @param jx splice junctions as GRanges objects
#'
#' @return a \code{\link[GenomicRanges]{GRangesList}} with altered transcripts.
#'
#' The algorithm works in the following steps.
#'
#' - Get transcripts for which the whole transcript range contains `pos1` and
#' `pos2` as GRangesList of GRanges of exons
#' - For each tx in transcripts do
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
#' Limitations:
#'  - [ ] vectorize over multiple input junctions
#'  - [ ] Both junction positions in same exon (Exitrons) are ignored
#'  - [ ] consider multiple junctions in same transcripts
#'
#' @examples
#'
#' @export
add_junc <- function(tx, jx){

  ## assume same length of tx and jx
  stopifnot(length(tx) == length(jx))
  stopifnot(!is.null(names(tx)))

  # get individual GRanges objects for start and end position of junction
  j1 <- GenomicRanges::resize(jx, width = 1, fix="start")
  j2 <- GenomicRanges::resize(jx, width = 1, fix="end")

  # convert GRangesLists as natural list objects
  tx_lst <- as.list(tx)
  jx_lst <- S4Vectors::split(jx, 1:length(jx)) %>% as.list()
  j1_lst <- S4Vectors::split(j1, 1:length(j1)) %>% as.list()
  j2_lst <- S4Vectors::split(j2, 1:length(j2)) %>% as.list()
  # tx_id <- names(tx)

  # Get indexes of affected exons ==============================================
  # get exon index in transcripts that are affected by the first position of the junction
  exon_idx1 = purrr::map2_int(tx_lst, j1_lst, ~GenomicRanges::findOverlaps(.y, .x, select = "first"))

  # get index of next exon left of first junction position
  left_exon_idx = purrr::map2_int(tx_lst, j1_lst, ~GenomicRanges::follow(.y, .x, ignore.strand = TRUE))

  # get exon index in transcripts that are affected by the second position of the junction
  exon_idx2 = purrr::map2_int(tx_lst, j2_lst, ~GenomicRanges::findOverlaps(.y, .x, select = "first"))

  # get index of next exon right of second junction position
  right_exon_idx = purrr::map2_int(tx_lst, j2_lst, ~GenomicRanges::precede(.y, .x, ignore.strand = TRUE))

  # get the indexes of exons that needs to be modified.
  # If position is not on an exon, the index is NA and the left (for junction
  # pos 1) or right (for junction pos 2) exon needs to be considered.
  used_exon_idx1 = ifelse(!is.na(exon_idx1), exon_idx1, left_exon_idx)
  used_exon_idx2 = ifelse(!is.na(exon_idx2), exon_idx2, right_exon_idx)

  # Modify affected exon start and end coordinates =============================
  # modify exon end position for first junction pos
  tx_lst = grl_update_end_at(tx_lst, used_exon_idx1, BiocGenerics::start(jx))

  # modify exon start position for second junction pos
  tx_lst = grl_update_start_at(tx_lst, used_exon_idx2, BiocGenerics::end(jx))

  ## OLD code:
  # modify exon end position for first junction pos
  tx_lst = purrr::modify2(tx_lst, used_exon_idx1, ~update_end_at(.x, BiocGenerics::start(jx), .y))

  # # modify exon start position for second junction pos
  tx = purrr::modify2(tx, used_exon_idx2, ~update_start_at(.x, pos2, .y))

  # delete exon_rank colum
  tx = purrr::map(tx, function(tx){
    tx$exon_rank = NULL
    return(tx)
  })

#
#   # In case of pos1 and pos2 on the same exon, ignore transcript
#   tx_df <- tx_df %>%
#     # either different exon, or one position not on any exon
#     filter(exon_idx1 != exon_idx2 | is.na(exon_idx1) | is.na(exon_idx2), )

  # Modify affected exon start and end coordinates ============================
  tx_df <- tx_df %>%
    mutate(

      # modify exon end position for first junction pos
      tx = purrr::modify2(tx, used_exon_idx1, ~update_end_at(.x, pos1, .y)),

      # # modify exon start position for second junction pos
      tx = purrr::modify2(tx, used_exon_idx2, ~update_start_at(.x, pos2, .y)),

      # delete exon_rank colum
      tx = purrr::map(tx, function(tx){
        tx$exon_rank = NULL
        return(tx)
      }),
    )
  # remove all exons inbetween =================================================
  tx_df <- tx_df %>%
    mutate(
      tx = purrr::map(tx, ~.x[!IRanges::overlapsAny(.x, junc_space_gr, type="within")])
    )

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

  purrr::pmap(
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

  purrr::pmap(
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
