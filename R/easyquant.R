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
#' @return A tibble with with the re-quantification results. This tibble has the columns
#' - `name`: name of the input sequence.
#' -  `interval_position`: comma-separated start end position of the interval relative to input sequence.
#' - `overlap_interval_end_reads`: reads overlapping the the interval end.
#' - `span_interval_end_reads`: read pairs spanning the interval end.
#' - `within_interval_reads`: number of reads that map into the interval.
#' - `interval_coverage_perc`: interval coverage defined as the percentage of the given interval that is covered by reads.
#' - `interval_coverage_mean`: interval coverage defined as the mean number of reads covering a position in the given interval.
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
    readr::read_delim(delim = "\t") %>%
  return(dat_easyqant)
}


#' Maps the re-quantification result from EasyQuant. on the junction-transcript centric tibble by hash id.
#'
#' @param path_to_easyquant_folder The path to EasyQuant folder
#' @param junc_tib The junction-transcript centric tibble, i.e. each row represents an altered transcript. Must contain a column `hash_id` with hash ids that relate to the column `name` in the Easyquant output.
#'
#' @return Extended junction-transcript tibble with re-quantification results.  The following columns are added:
#' -  `interval_position`: comma-separated start end position of the interval relative to input sequence.
#' - `overlap_interval_end_reads`: reads overlapping the the interval end.
#' - `span_interval_end_reads`: read pairs spanning the interval end.
#' - `within_interval_reads`: number of reads that map into the interval.
#' - `interval_coverage_perc`: interval coverage defined as the percentage of the given interval that is covered by reads.
#' - `interval_coverage_mean`: interval coverage defined as the mean number of reads covering a position in the given interval.
#
#'
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


