
#' convert a `junc_id` into a dataset of altered transcript
#' \code{\link[GenomicRanges]{GRangesList}}.
#'
#' @param chr chromosome
#' @param pos1 start position of junction
#' @param pos2 start position of junction
#' @param transcripts a GRangesList with transcripts defined as GRanges of exons
#'   created by `GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)`.
#'
#' @return a GRangesList with alternative transcripts.
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
#'   - remove all exons completeley contained in the interval [`pos1`, `pos2`].
#'
#' TODO:
#'  - [ ] vectorize over multiple input junctions
#'  - [x] add `tx_name` as output column
#'  - [ ] Support junctions with both positions in the same exon
#'
#' @examples
#'
#' junc_to_tx("chr2", 152389996, 152392205, toy_transcripts)
#'
#' @export
junc_to_tx <- function(chr, pos1, pos2, transcripts){

  ## assume pos1 to be smaller than pos2
  stopifnot(length(chr) == 1)
  stopifnot(length(pos1) == 1)
  stopifnot(length(pos2) == 1)
  # stopifnot(length(strand) == 1)
  stopifnot(pos1 < pos2)

  gr1 <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos1, pos1))
  gr2 <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos2, pos2))

  junc_space_gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(pos1+1, pos2-1))

  # get subset of transcripts that overlap with both positions =================
  tx_ranges <- base::range(transcripts)

  tx_idx1 <- GenomicRanges::findOverlaps(gr1, tx_ranges) %>% S4Vectors::subjectHits()
  tx_idx2 <- GenomicRanges::findOverlaps(gr2, tx_ranges) %>% S4Vectors::subjectHits()
  txs <- transcripts[intersect(tx_idx1, tx_idx2)]

  # Get indexes of affected exons ==============================================
  tx_df <- tibble::tibble(
    tx = as.list(txs),
    tx_name = names(txs),
    exon_idx1 = purrr::map_int(tx, ~GenomicRanges::findOverlaps(gr1, .x, select = "first")),
    left_exon_idx = purrr::map_int(tx, ~GenomicRanges::follow(gr1, .x, ignore.strand = TRUE)),
    exon_idx2 = purrr::map_int(tx, ~GenomicRanges::findOverlaps(gr2, .x, select = "first")),
    right_exon_idx = purrr::map_int(tx, ~GenomicRanges::precede(gr2, .x, ignore.strand = TRUE)),

    # get the indexes of exons that needs to be modified
    used_exon_idx1 = ifelse(!is.na(exon_idx1), exon_idx1, left_exon_idx),
    used_exon_idx2 = ifelse(!is.na(exon_idx2), exon_idx2, right_exon_idx),
  )

  # In case of pos1 and pos2 on the same exon, ignore transcript
  tx_df <- tx_df %>%
    # either different exon, or one position not on any exon
    filter(exon_idx1 != exon_idx2 | is.na(exon_idx1) | is.na(exon_idx2), )

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

  # add coordinate in transcript sequence ====================================
  tx_grl <- GenomicRanges::GRangesList(tx_df$tx)
  names(tx_grl) <- tx_df$tx_name
  pos1_tx <- GenomicFeatures::mapToTranscripts(gr1, tx_grl) %>% BiocGenerics::start()
  pos2_tx <- GenomicFeatures::mapToTranscripts(gr2, tx_grl) %>% BiocGenerics::start()
  strand_tx <- tx_grl %>% BiocGenerics::strand() %>% unique() %>% as.character()
  junc_pos_tx <- ifelse(strand_tx == "+", pos1_tx, pos2_tx)

  tx_df <- tx_df %>%
    mutate(
      junc_pos_tx = junc_pos_tx
    )

  return(tx_df)
}


#' update end coordinate of a subset of ranges by given indices
#'
#' @param gr A GRanges object
#' @param pos The new coordinate
#' @param at The indices of ranges to update
#'
#' @return A modified GRanges object
update_end_at <- function(gr, pos, at){
  BiocGenerics::end(gr[at]) <- pos
  return(gr)
}

#' update start coordinate of a subset of ranges by given indices
#'
#' @param gr A GRanges objec
#' @param pos The new coordinate
#' @param at The indeces of ranges to update
#'
#' @return A modified GRanges object
update_start_at <- function(gr, pos, at){
  BiocGenerics::start(gr[at]) <- pos
  return(gr)
}
