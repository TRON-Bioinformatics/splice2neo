#' Creates a table with context sequences for re-quantification analysis with EasyQuant.
#'
#' @param df  A data.frame with splice junctions in rows and at least the columns:
#'
#'   -  `cts_id` hash id for a context sequence
#'   -  `cts_seq` The transcript context sequence.
#'   -  `cts_junc_pos` The position of the junction within the context sequence
#'
#'
#' @return A tibble that can be used for re-quantification with EasyQuant This tibble has the columns
#' - `name`: Name of the input sequence.
#' -  `sequence`: The context sequence that will be requantified.
#' - `position`: The position of the junction within the context sequence.
#'
#' @examples
#'
#'
#'@import dplyr
#'
#'@export
transform_for_requant <- function(df){

  stopifnot("cts_id" %in% names(df))
  stopifnot("cts_seq" %in% names(df))
  stopifnot("cts_junc_pos" %in% names(df))

  dat <- df %>%
    dplyr::select(cts_id, cts_seq, cts_junc_pos) %>%
    dplyr::rename(name = cts_id,
                  sequence = cts_seq,
                  position = cts_junc_pos)
  dat <- dat %>%
    distinct()%>%
    filter(!is.na(sequence))
  return(dat)
}


#' Imports the re-quantification results from analysis with EasyQuant
#'
#' @param path_folder The path to the EasyQuant folder
#'
#' @return A tibble with with the re-quantification results. More details about Easyquant can be obtained at https://github.com/TRON-Bioinformatics/easyquant
#'
#' @examples
#'
#'
#'@import dplyr
#'
#'@export
read_requant <- function(path_folder){
  path.to.easyquant.file <- paste(path_folder, "quantification.tsv" ,sep = "/" )
  if(!file.exists(path.to.easyquant.file)){
    stop("quantification.tsv file is missing")
  }
  dat_easyqant <- path.to.easyquant.file %>%
    readr::read_delim(delim = "\t")

  # intron retentions
  # interval-1: end --> exon-intron boundary
  # interval-2: interval --> intron; end --> exon-intron boundary
  dat_ir <- dat_easyqant %>%
    group_by(name) %>%
    filter(n() == 3)

  if(nrow(dat_ir) > 0){
    dat_ir <-  dat_ir %>%
      mutate(interv = c(1, 2, 3)) %>%
      pivot_wider(
        names_from = interv,
        values_from = c(
          overlap_interval_end_reads,
          span_interval_end_pairs,
          within_interval,
          coverage_perc,
          coverage_mean,
          coverage_median,
          interval
        )
      )

    dat_ir <- dat_ir %>%
      dplyr::select(-overlap_interval_end_reads_3, -span_interval_end_pairs_3,) %>%
      dplyr::rename(
        junc_interval_start = overlap_interval_end_reads_1,
        junc_interval_end = overlap_interval_end_reads_2,
        span_interval_start = span_interval_end_pairs_1,
        span_interval_end = span_interval_end_pairs_2,
        within_interval = within_interval_2,
        within_interval_left = within_interval_1,
        within_interval_right = within_interval_3,
        coverage_perc = coverage_perc_2,
        coverage_perc_left = coverage_perc_1,
        coverage_perc_right = coverage_perc_3,
        coverage_mean = coverage_mean_2,
        coverage_mean_left = coverage_mean_1,
        coverage_mean_right = coverage_mean_3,
        coverage_median = coverage_median_2,
        coverage_median_left = coverage_median_1,
        coverage_median_right = coverage_median_3,
        interval = interval_2,
        interval_left = interval_1,
        interval_right = interval_3
      )
  } else {
    dat_ir <- dat_ir %>%
      mutate(junc_interval_end = NA,
             span_interval_end = NA,
             within_interval = NA,
             coverage_median = NA,
             interval = NA)
  }


  # non- intron retentions
  # interval-1: end --> exon-intron boundary / junction of interest
  # add interval info at interval_left / interval_right
  dat_no_ir <- dat_easyqant %>%
    group_by(name) %>%
    filter(n() == 2)

  if(nrow(dat_no_ir) > 0){
    dat_no_ir <- dat_no_ir%>%
      mutate(interv = c(1, 2))%>%
      pivot_wider(
        names_from = interv,
        values_from = c(
          overlap_interval_end_reads,
          span_interval_end_pairs,
          within_interval,
          coverage_perc,
          coverage_mean,
          coverage_median,
          interval
        )
      )
    dat_no_ir <- dat_no_ir %>%
      dplyr::rename(
        junc_interval_start = overlap_interval_end_reads_1,
        span_interval_start = span_interval_end_pairs_1,
        within_interval_left = within_interval_1,
        within_interval_right = within_interval_2,
        coverage_perc_left = coverage_perc_1,
        coverage_perc_right = coverage_perc_2,
        coverage_mean_left = coverage_mean_1,
        coverage_mean_right = coverage_mean_2,
        coverage_median_left = coverage_median_1,
        coverage_median_right = coverage_median_2,
        interval_left = interval_1,
        interval_right = interval_2
      ) %>%
      select(-overlap_interval_end_reads_2, -span_interval_end_pairs_2 )
  }

  dat_easyqant <- bind_rows(dat_ir, dat_no_ir) %>%
    ungroup() %>%
    select(name,
           junc_interval_start, junc_interval_end,
           span_interval_start, span_interval_end,
           within_interval, within_interval_left, within_interval_right,
           coverage_perc, coverage_perc_left, coverage_perc_right,
           coverage_mean, coverage_mean_left, coverage_mean_right,
           coverage_median, coverage_median_left, coverage_median_right,
           interval, interval_left, interval_right)

  return(dat_easyqant)
}


#' Maps the re-quantification result from EasyQuant. on the junction-transcript centric tibble by hash id.
#'
#' @param path_to_easyquant_folder The path to EasyQuant folder
#' @param junc_tib The junction-transcript centric tibble, i.e. each row represents an altered transcript. Must contain a column `hash_id` with hash ids that relate to the column `name` in the Easyquant output.
#'
#' @return Extended junction-transcript tibble with re-quantification results by Easyquant. More details about Easyquant can be obtained at https://github.com/TRON-Bioinformatics/easyquant
#' The following columns are added:
#' - `junc_interval_start`: Junction reads overlapping the junction of interest.
#'  In case of intron retentions, junctions reads that overlap the start of the intron retention interval.
#' - `junc_interval_end`: In case of intron retentions, junctions reads that overlap the end of the intron retention interval.
#' Is NA for non intron retention events.
#' - `span_interval_start`: Spanning pairs overlapping the junction of interest. In case of intron retentions, spanning pairs that overlap the start of the intron retention interval.
#' - `span_interval_end`: In case of intron retentions, spanning pairs that overlap the end of the intron retention interval.
#' Is NA for non intron retention events.
#' - `within_interval`: Number of reads that map into the intron retention interval.
#' Is NA for non intron retention events.
#' - `within_interval_left`: Number of reads that map into interval left to the intron retention or splice junction of interest.
#' - `within_interval_right`: Number of reads that map into interval right of the intron retention or splice junction of interest.
#' - `coverage_perc`: Intron retention interval coverage defined as the percentage of the given interval that is covered by reads.
#' Is NA for non intron retention events.
#' - `coverage_perc_left`: Interval coverage defined as the percentage of the interval
#' left to the intron retention or splice junction of interest that is covered by reads.
#' - `coverage_perc_right`: Interval coverage defined as the percentage of the interval
#' right to the intron retention or splice junction of interest that is covered by reads.
#' - `coverage_mean`: Interval coverage defined as the mean number of reads covering a position in the given interval representing an intron retention.
#' Is NA for non intron retention events.
#' - `coverage_mean_left`: Interval coverage defined as the mean number of reads covering a position in the interval left to the intron retention or splice junction of interest.
#' - `coverage_mean_right`: Interval coverage defined as the mean number of reads covering a position in the interval right to the intron retention or splice junction of interest.
#' - `coverage_median`: Interval coverage defined as the median number of reads covering a position in the given interval representing an intron retention. Is NA for non intron retention events.
#' - `coverage_median_left`: Interval coverage defined as the median number of reads covering a position in the interval left to the intron retention or splice junction of interest.
#' - `coverage_median_right`: Interval interval coverage defined as the median number of reads covering a position in the interval right to the intron retention or splice junction of interest.
#' @examples
#'
#'
#'@import dplyr
#'
#'@export
map_requant <- function(path_to_easyquant_folder, junc_tib) {

  stopifnot("cts_id" %in% names(junc_tib))

  dat_easyqant <- read_requant(path_to_easyquant_folder)
  dat_junc <- junc_tib %>%
    left_join(dat_easyqant, by = c("cts_id" = "name"))
  return(dat_junc)
}


