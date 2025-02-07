#' Annotate if splice junction points overlap with genomic positions of retroelements.
#'
#' @param df A data.frame with splice junctions in rows and at least the columns:
#'
#'   -  `junc_id` junction id consisting of genomic coordinates
#'
#' @param rmsk \code{\link[GenomicRanges]{GRanges}} of retroelements (RepeatMasker)
#'
#' @return A data.frame as the input with additional columns annotating overlaps with retroelemts.
#'
#'  - `potential_jet` Indicator if junction is potentially a jet
#'  - `pos_left_retroelement`: Name of retroelement overlapping left splice site
#'  - `pos_right_retroelement`: Name of retroelement overlapping right splice site
#'
#' @examples
#' ah <- AnnotationHub()
#' query(ah, c("RepeatMasker", "Homo sapiens"))
#' rmsk <- ah[["AH99003"]]
#'
#' junc_df <- tibble::tibble(
#'   junc_id = c("chr2:152389996-152392205:-", "chr2:152389996-152390729:-",
#'               "chr2:152389955-152389956:-")
#' )
#'
#' annotate_potential_jet(junc_df, rmsk)
#'
#' @import dplyr
#' @export
annotate_potential_jet <- function(df, rmsk) {
  
  stopifnot(is.data.frame(df))
  stopifnot("junc_id" %in% names(df))
  
  # build GRanges of junctions
  jx <- junc_to_gr(df$junc_id)
  
  # junction starts
  junc_starts <- GenomicRanges::resize(jx, width = 1, fix = "start")
  junc_starts_grl <- S4Vectors::split(junc_starts, 1:length(jx))
  
  # junction ends
  junc_ends <- GenomicRanges::resize(jx, width = 1, fix = "end")
  junc_ends_grl <- S4Vectors::split(junc_ends, 1:length(jx))
  
  df <- df %>%
    dplyr::mutate(
        pos_left_lst  = junc_starts_grl %>% as.list(),
        pos_right_lst = junc_ends_grl %>% as.list()
    )
  
  df <- df %>%
    dplyr::mutate(
      pos_left_retroelement_idx =
        purrr::map_dbl(pos_left_lst, ~ GenomicRanges::findOverlaps(.x, rmsk, select = "first")),
      pos_right_retroelement_idx =
        purrr::map_dbl(pos_right_lst, ~ GenomicRanges::findOverlaps(.x, rmsk, select = "first")),
      pos_left_retroelement =
        dplyr::if_else(is.na(pos_left_on_retroelement), NA, rmsk$repName[pos_left_retroelement_idx]),
      pos_right_retroelement =
        dplyr::if_else(is.na(pos_right_on_retroelement), NA, rmsk$repName[pos_right_retroelement_idx]),
    ) %>%
    mutate(
      potential_jet = dplyr::if_any(c(pos_left_retroelement, pos_right_retroelement),~!is.na(.x))
    )
  
  df <- df %>%
    dplyr::select(-pos_left_lst, -pos_right_lst, -pos_left_retroelement_idx, -pos_right_retroelement_idx)
  
  return(df)

}
