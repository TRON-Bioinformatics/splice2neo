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

  dat_easyqant <- dat_easyqant %>%
    group_by(name) %>%
    filter(if(n() == 2) row_number() == 1 else row_number() != 3)

  # intron retentions
  # interval-1: end --> exon-intron boundary
  # interval-2: interval --> intron; end --> exon-intron boundary
  dat_ir <- dat_easyqant %>%
    filter(n() > 1) %>%
    mutate(interv = c(1,2)) %>%
    pivot_wider(names_from = interv, values_from = c(overlap_stop, span_read, within_interval, coverage_perc, coverage_mean, interval))

  # non- intron retentions
  # interval-1: end --> exon-intron boundary / junction of interest
  # within_interval/coverage not relevant for these events --> set to NA
  dat_no_ir <- dat_easyqant %>%
    filter(n() == 1) %>%
    select(-interval) %>%
    dplyr::rename(junc_interval_start = overlap_stop) %>%
    dplyr::rename(span_interval_start = span_read) %>%
    mutate(within_interval = NA) %>%
    mutate(coverage_perc = NA)%>%
    mutate(coverage_mean = NA)


  dat_ir <- dat_ir %>%
    dplyr::select(-within_interval_1, -coverage_perc_1,-coverage_mean_1) %>%
    dplyr::rename("junc_interval_start" = "overlap_stop_1") %>%
    dplyr::rename("junc_interval_end" = "overlap_stop_2") %>%
    dplyr::rename("span_interval_start" = "span_read_1") %>%
    dplyr::rename("span_interval_end" = "span_read_2") %>%
    dplyr::rename("within_interval" = "within_interval_2") %>%
    dplyr::rename("coverage_perc" = "coverage_perc_2")%>%
    dplyr::rename("coverage_mean" = "coverage_mean_2") %>%
    dplyr::select(-interval_1, - interval_2)

  dat_easyqant <- bind_rows(dat_ir, dat_no_ir)

  return(dat_easyqant)
}


#' Maps the re-quantification result from EasyQuant. on the junction-transcript centric tibble by hash id.
#'
#' @param path_to_easyquant_folder The path to EasyQuant folder
#' @param junc_tib The junction-transcript centric tibble, i.e. each row represents an altered transcript. Must contain a column `hash_id` with hash ids that relate to the column `name` in the Easyquant output.
#'
#' @return Extended junction-transcript tibble with re-quantification results by Easyquant. More details about Easyquant can be obtained at https://github.com/TRON-Bioinformatics/easyquant
#' The following columns are added:
#' - `junc_interval_start`: Junction reads overlapping the junction of interest. In case of intron retentions, junctions reads that overlap the start of the intron retention interval.
#' - `junc_interval_end`: In case of intron retentions, junctions reads that overlap the end of the intron retention interval. Is NA for non intron retention events.
#' - `span_interval_start`: Spanning pairs overlapping the junction of interest. In case of intron retentions, spanning pairs that overlap the start of the intron retention interval.
#' - `span_interval_end`: In case of intron retentions, spanning pairs that overlap the end of the intron retention interval. Is NA for non intron retention events.
#' - `within_interval`: Number of reads that map into the intron retention interval. Is NA for non intron retention events.
#' - `coverage_perc`: Intron retention interval coverage defined as the percentage of the given interval that is covered by reads. Is NA for non intron retention events.
#' - `coverage_mean`: Intron retention interval coverage defined as the mean number of reads covering a position in the given interval. Is NA for non intron retention events.
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


