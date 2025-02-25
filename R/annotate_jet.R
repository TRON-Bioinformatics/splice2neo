#' Annotate if splice sites of a junction overlap with genomic positions of transposable elements.
#'
#' @param df A data.frame with splice junctions in rows and at least the column:
#'
#'   -  `junc_id` junction id consisting of genomic coordinates
#'
#' @param rmsk \code{\link[GenomicRanges]{GRanges}} of transposable elements (e.g. RepeatMasker)
#'
#' @return A data.frame as the input with additional columns annotating overlaps with transposable elements.
#'
#'  - `potential_jet` Indicator if junction is potentially a JET
#'  - `left_side_retroelement`: Name of transposable element overlapping left splice site
#'  - `left_side_retroelement_class`: Class of transposable element overlapping left splice site
#'  - `right_side_retroelement`: Name of transposable element overlapping right splice site
#'  - `right_side_retroelement_class`: Class of transposable element overlapping right splice site
#'
#' @examples
#'
#'
#' rmsk <- readr::read_tsv(system.file("extdata", "rmsk_hg19_subset.tsv.gz", package = "splice2neo"))
#' rmsk <- GenomicRanges::makeGRangesFromDataFrame(rmsk, keep.extra.columns = T)
#' junc_df <- tibble::tibble(
#'   junc_id = c("chr2:152407458-152408252:-")
#' )
#'
#' annotate_potential_jet(junc_df, rmsk)
#'
#' @import dplyr purrr GenomicRanges tibble
#' @export
annotate_potential_jet <- function(df, rmsk) {

  stopifnot(is.data.frame(df))
  stopifnot(class(rmsk) %in% c("GRanges", "CompressedGRanges"))
  stopifnot("junc_id" %in% names(df))

  # build GRanges of junctions
  jx <- junc_to_gr(df$junc_id)

  # junction starts
  junc_starts <- GenomicRanges::resize(jx, width = 1, fix = "start", ignore.strand = TRUE)

  # junction ends
  junc_ends <- GenomicRanges::resize(jx, width = 1, fix = "end", ignore.strand = TRUE)

  # Annotate overlaps of left splice-site with RepeatMasker predictions
  left_overlaps <- GenomicRanges::findOverlaps(junc_starts, rmsk, ignore.strand=F)
  left_jx_idx <- left_overlaps %>% S4Vectors::queryHits()
  left_retroelement_idx <- left_overlaps %>% S4Vectors::subjectHits()

  left_to_retroelement <- tibble::tibble(
    junc_id = df$junc_id[left_jx_idx],
    pos_left_retroelement = rmsk$repName[left_retroelement_idx],
    pos_left_retroelement_class = rmsk$repClass[left_retroelement_idx]
  )

  # Annotate overlaps of right splice-site with RepeatMasker predictions
  right_overlaps <- GenomicRanges::findOverlaps(junc_ends, rmsk, ignore.strand=F)
  right_jx_idx <- right_overlaps %>% S4Vectors::queryHits()
  right_retroelement_idx <- right_overlaps %>% S4Vectors::subjectHits()

  right_to_retroelement <- tibble::tibble(
    junc_id = df$junc_id[right_jx_idx],
    pos_right_retroelement = rmsk$repName[right_retroelement_idx],
    pos_right_retroelement_class = rmsk$repClass[right_retroelement_idx]
  )

  # Combine by junc_id. If multiple Elements overlap they are collapsed
  retro_element_overlap <- left_to_retroelement %>%
    full_join(right_to_retroelement) %>%
    dplyr::distinct() %>%
    dplyr::group_by(junc_id) %>%
    dplyr::summarize(
      left_side_retroelement = str_c(unique(purrr::discard(pos_left_retroelement, is.na)), collapse = "|"),
      left_side_retroelement_class = str_c(unique(purrr::discard(pos_left_retroelement_class, is.na)), collapse = "|"),
      right_side_retroelement = str_c(unique(purrr::discard(pos_right_retroelement, is.na)), collapse = "|"),
      right_side_retroelement_class = str_c(unique(purrr::discard(pos_right_retroelement_class, is.na)), collapse = "|"),
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate_all(., list(~dplyr::na_if(.,"")))

  df <- df %>%
    dplyr::left_join(retro_element_overlap) %>%
    dplyr::mutate(
      potential_jet = dplyr::if_any(c(left_side_retroelement, right_side_retroelement), ~!is.na(.x))
    )

  return(df)
}
