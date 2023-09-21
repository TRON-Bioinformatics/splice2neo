
#' Select transcripts per junction that are more likely to be affected by a junction.
#'
#' @param df A data.frame with splice junctions in rows and at least the columns:
#'
#'   -  `junc_id` junction id consisting of genomic coordinates
#'   -  `tx_id` junction id consisting of genomic coordinates
#'   -  `tx_lst` a list of \code{\link[GenomicRanges]{GRanges}} with the transcript
#'
#'
#' @return A data.frame as the input but with relevant transcript and junction combinations
#' NOTE: This function selects transcripts that are more likely to be affected to reduce data dimension.
#' It excludes transcripts in which both junction positions are located in an intron. Junctions in a given transcript must be either be classified as
#' ES, IR, exitron or ASS event or have both junction positions in an exon. Other junction-transcript combinations are also excluded.
#' -
#'
#'
#' @import dplyr
#' @export
choose_tx <- function(df){

  stopifnot(is.data.frame(df))
  stopifnot("junc_id" %in% names(df))
  stopifnot("tx_id" %in% names(df))
  stopifnot("tx_lst" %in% names(df))

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

  df <- df %>%
    tidyr::separate(junc_id, into = c("chr","start-end","strand" ), sep = ":", remove = FALSE) %>%
    tidyr::separate(`start-end`, into = c("start", "end" ), sep = "-") %>%
    mutate(
      start = as.numeric(start),
      end = as.numeric(end)
    ) %>%
    mutate(
      exon_starts_lst  = as.list(exon_starts),
      exon_ends_lst = as.list(exon_ends),
      pos_left_lst  = junc_starts_grl %>% as.list() ,
      pos_right_lst = junc_ends_grl %>% as.list()

    )

  # annotate the exon details for each junction - transcript combination
  df_annot <- annotate_exon_idx(df)
  # annotate putative classes of junction - transcripts
  df_annot_class <- classify_transcripts(df_annot)

  df_fil <- df_annot_class %>%
    filter(!is.na(putative_event_type))


  # join the original input data.frame with the association of junction to transcripts
  out_df <- df %>%
    dplyr::left_join(junc_to_tx, by = c("junc_id", "tx_id"))

  return(out_df)
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
    mutate(intron_retention = ifelse(abs(start - end) == 1, TRUE, FALSE)) %>%
    mutate(pos_left_on_exon = ifelse(!is.na(pos_left_on_exon_idx) , TRUE, FALSE)) %>%
    mutate(pos_right_on_exon = ifelse(!is.na(pos_right_on_exon_idx) , TRUE, FALSE)) %>%
    mutate(pos_left_on_exon_boundary = ifelse(!is.na(pos_left_on_exon_end) , TRUE, FALSE)) %>%
    mutate(pos_right_on_exon_boundary = ifelse(!is.na(pos_right_on_exon_start) , TRUE, FALSE))



  df2 <- df1 %>%
    group_by(junc_id) %>%
    mutate(putative_event_type =
             case_when(
               # ES
               any(both_pos_on_exon_boundaries) & both_pos_on_exon_boundaries ~ "ES",
               any(both_pos_on_exon_boundaries) & !both_pos_on_exon_boundaries ~ "rm",
               # do not use junctions with both coordinates in introns
               # TODO: discuss if relevant
               !pos_left_on_exon & !pos_right_on_exon ~ "rm",
               # exitrons
               both_pos_not_on_exon_boundaries & both_pos_on_same_exon ~ "exitron",
               both_pos_not_on_exon_boundaries & any_pos_on_intron ~ "rm",
               both_pos_not_on_exon_boundaries ~ "complex",
               # IR on the right side
               pos_left_on_exon_boundary & intron_retention ~ "IR",
               # ASS on the right side
               pos_left_on_exon_boundary ~ "ASS_right",
               # IR on the left side
               pos_right_on_exon_boundary & intron_retention ~ "IR",
               # ASS on the left side
               pos_right_on_exon_boundary ~ "ASS_left",
               TRUE ~ "rm"
               # TODO: how to generalise for jets
             )) %>%
    mutate(putative_event_type = ifelse(putative_event_type == "rm", NA, putative_event_type))

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
          as.list(tx_glst),
          ~ GenomicRanges::findOverlaps(.x, .y, select = "first")
        ),
      # junction start end idx
      pos_right_on_exon_idx =
        purrr::map2_dbl(
          pos_right_lst,
          as.list(tx_glst),
          ~ GenomicRanges::findOverlaps(.x, .y, select = "first")
        )
    )

}

