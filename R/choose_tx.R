
#' Select transcripts per junction that are more likely to be affected by a junction. All possible transcripts can be annotated with `add_tx()`
#'
#' @param df A data.frame with splice junctions in rows and at least the columns:
#'
#'   -  `junc_id` junction id consisting of genomic coordinates
#'   -  `tx_id` junction id consisting of genomic coordinates
#'   -  `tx_lst` a list of \code{\link[GenomicRanges]{GRanges}} with the transcript
#'
#'
#' @return A data.frame as with relevant transcript and junction combinations. If tx_id is NA in the input data.frame such rows are removed from the output data.frame.
#'
#'
#' NOTE: This function selects transcripts that are more likely to be affected to reduce data dimension.
#' It excludes transcripts in which both junction positions are located in an intron. Junctions in a given transcript must be either be classified as
#' ES, IR, exitron or ASS event or have both junction positions in an exon. Other junction-transcript combinations are also excluded.
#' This function may loose relevant or keep irrelevant junction-transcripts in particular in regions with mutliple isoforms with distinct splicing pattern.
#'
#'
#'
#' @import dplyr
#' @export
choose_tx <- function(df){

  stopifnot(is.data.frame(df))
  stopifnot("junc_id" %in% names(df))
  stopifnot("tx_id" %in% names(df))
  stopifnot("tx_lst" %in% names(df))

  df <- df %>%
    filter(!is.na(tx_id))

  # build GRanges of junction
  jx <- junc_to_gr(df$junc_id)

  # convert GRangesLists as natural list objects
  tx_glst <- GenomicRanges::GRangesList(df$tx_lst)

  # exon starts
  exon_starts <- GenomicRanges::resize(tx_glst, width = 1, fix="start")

  # exon ends
  exon_ends <- GenomicRanges::resize(tx_glst, width = 1, fix="end")

  # junction starts
  junc_starts <- GenomicRanges::resize(jx, width = 1, fix="start")
  junc_starts_grl <- S4Vectors::split(junc_starts, 1:length(jx))

  # junction ends
  junc_ends <- GenomicRanges::resize(jx, width = 1, fix="end")
  junc_ends_grl <- S4Vectors::split(junc_ends, 1:length(jx))

  df1 <- df %>%
    tidyr::separate(junc_id, into = c("chr_junc","start-end","strand_junc" ), sep = ":", remove = FALSE) %>%
    tidyr::separate(`start-end`, into = c("start_junc", "end_junc" ), sep = "-") %>%
    mutate(
      start_junc = as.numeric(start_junc),
      end_junc = as.numeric(end_junc)
    ) %>%
    mutate(
      exon_starts_lst  = as.list(exon_starts),
      exon_ends_lst = as.list(exon_ends),
      pos_left_lst  = junc_starts_grl %>% as.list() ,
      pos_right_lst = junc_ends_grl %>% as.list()

    )

  # annotate the exon details for each junction - transcript combination
  df_annot <- annotate_exon_idx(df1)
  # annotate putative classes of junction - transcripts
  df_annot_class <- classify_transcripts(df_annot)

  df_fil <- df_annot_class %>%
    filter(!is.na(putative_event_type)) %>%
    dplyr::select(
      - chr_junc,
      - start_junc,
      - end_junc,
      - strand_junc,
      - exon_starts_lst,
      - exon_ends_lst,
      - pos_left_lst,
      - pos_right_lst,
      - pos_left_on_exon_end,
      - pos_right_on_exon_start,
      - pos_left_on_exon_idx,
      - pos_right_on_exon_idx,
      - both_pos_on_exon_boundaries,
      - both_pos_not_on_exon_boundaries,
      - both_pos_on_same_exon,
      - any_pos_on_intron,
      - intron_retention,
      - pos_left_on_exon,
      - pos_right_on_exon,
      - pos_left_on_exon_boundary,
      - pos_right_on_exon_boundary
    )

  return(df_fil)
}


#' Classify splice junctions- transcripts combinations based on junction and coordinates of annotated exons
#'
#' @param df A data.frame with unctions- transcripts combination in rows.
#'#'
#'
#' @return A data.frame as the input but with more annotations.
#'
#'#'
#' @import dplyr
classify_transcripts <- function(df){

  df1 <- df %>%
    mutate(both_pos_on_exon_boundaries = ifelse(
      !is.na(pos_left_on_exon_end) &
        !is.na(pos_right_on_exon_start),
      TRUE,
      FALSE
    )) %>%
    mutate(both_pos_not_on_exon_boundaries = ifelse(
      is.na(pos_left_on_exon_end) &
        is.na(pos_right_on_exon_start),
      TRUE,
      FALSE
    )) %>%
    mutate(
      both_pos_on_same_exon = ifelse(pos_left_on_exon_idx == pos_right_on_exon_idx, TRUE, FALSE)
    ) %>%
    mutate(any_pos_on_intron = ifelse(
      is.na(pos_left_on_exon_idx) |
        is.na(pos_right_on_exon_idx),
      TRUE,
      FALSE
    )) %>%
    mutate(intron_retention = ifelse(abs(start_junc - end_junc) == 1, TRUE, FALSE)) %>%
    mutate(pos_left_on_exon = ifelse(!is.na(pos_left_on_exon_idx) , TRUE, FALSE)) %>%
    mutate(pos_right_on_exon = ifelse(!is.na(pos_right_on_exon_idx) , TRUE, FALSE)) %>%
    mutate(pos_left_on_exon_boundary = ifelse(!is.na(pos_left_on_exon_end) , TRUE, FALSE)) %>%
    mutate(pos_right_on_exon_boundary = ifelse(!is.na(pos_right_on_exon_start) , TRUE, FALSE))



  df2 <- df1 %>%
    group_by(junc_id) %>%
    mutate(putative_event_type =
             case_when(
               both_pos_on_exon_boundaries & abs(pos_left_on_exon_end - pos_right_on_exon_start) == 1 ~ "WT junction",
               # ES
               any(both_pos_on_exon_boundaries) & both_pos_on_exon_boundaries ~ "ES",
               # there are cases for which ES in one transcript can be an exitron in another transcript
               any(both_pos_on_exon_boundaries) & both_pos_on_same_exon ~ "exitron",
               any(both_pos_on_exon_boundaries) & abs(pos_left_on_exon_idx - pos_right_on_exon_idx) == 1 & !is.na(pos_right_on_exon_idx) ~ "ASS",
               any(both_pos_on_exon_boundaries) & !both_pos_on_exon_boundaries ~ "rm",
               # do not use junctions with both coordinates in introns
               # TODO: discuss if relevant
               !pos_left_on_exon & !pos_right_on_exon ~ "rm",
               # exitrons
               # check that this is not an "intron retention" (pos1:(pos1+1))
               both_pos_not_on_exon_boundaries & both_pos_on_same_exon & intron_retention ~ "rm",
               both_pos_not_on_exon_boundaries & both_pos_on_same_exon  ~ "exitron",
               both_pos_not_on_exon_boundaries & any_pos_on_intron ~ "rm",
               both_pos_not_on_exon_boundaries ~ "complex",
               # IR on the right side
               pos_left_on_exon_boundary & intron_retention ~ "IR",
               # ASS on the right side
               pos_left_on_exon_boundary ~ "ASS",
               # IR on the left side
               pos_right_on_exon_boundary & intron_retention ~ "IR",
               # ASS on the left side
               pos_right_on_exon_boundary ~ "ASS",
               TRUE ~ "rm"

             )) %>%
    mutate(putative_event_type = ifelse(putative_event_type == "rm", NA, putative_event_type)) %>%
    ungroup()

}



#' Annotate splice junctions- transcripts combination and transcripts with exon details
#'
#' @param df  A data.frame with junctions- transcripts combination in rows.
#'
#'   -  `junc_id` junction id consisting of genomic coordinates
#'   -  `tx_id` junction id consisting of genomic coordinates
#'   -  `tx_lst` a list of \code{\link[GenomicRanges]{GRanges}} with the transcript
#'
#'
#' @return A data.frame as the input but with more annotations.
#'
#'
#'
#' @import dplyr
annotate_exon_idx <- function(df){

  df %>%
    mutate(
      # junction start on exon end?
      pos_left_on_exon_end =
        purrr::map2_dbl(
          pos_left_lst,
          exon_ends_lst,
          ~ GenomicRanges::findOverlaps(.x, .y, select = "first")
        ),
      # junction end on exon start?
      pos_right_on_exon_start =
        purrr::map2_dbl(
          pos_right_lst,
          exon_starts_lst,
          ~ GenomicRanges::findOverlaps(.x, .y, select = "first")
        ),
      # junction start exon idx
      pos_left_on_exon_idx =
        purrr::map2_dbl(
          pos_left_lst,
          tx_lst,
          ~ GenomicRanges::findOverlaps(.x, .y, select = "first")
        ),
      # junction start end idx
      pos_right_on_exon_idx =
        purrr::map2_dbl(
          pos_right_lst,
          tx_lst,
          ~ GenomicRanges::findOverlaps(.x, .y, select = "first")
        )
    )

}

