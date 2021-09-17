#' Create table for re-quantification analysis with easyquant
#'
#' @param cts_id An id for the cts sequence to be re-quantified. If the
#' context sequence was determined with the function `junct_to_cts` the associated hash id `cts_id` can be used.
#' @param cts_seq The context sequence.
#' @param junc_position The position of the junction within the context sequence.  If the
#' context sequence was determined with `junct_to_cts` the associated hash id `cts_junc_pos` can be used.
#'
#' @return A tibble that can be used for requantification with easyquant. This tibble has the columns
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
transform_for_requant <- function(cts_id, cts_seq, junc_position){
  dat <- tibble(
    name = cts_id,
    sequence = cts_seq,
    position = junc_position
  )
  dat <- dat %>%
    distinct()%>%
    filter(!is.na(sequence))
  return(dat)
}


#' Imports re-quantification from analysis with easyquant
#'
#' @param path_folder The path to easyquant folder
#'
#' @return A tibble with with the requantification results. This tibble has the columns
#' - `name`: name of the input sequence. this will be a hash id.
#' -  `pos`: position of interest relative to input sequence
#' - `junc`: reads overlapping the position of interest
#' - `span`: read pairs spanning the position of interest
#' - `anch`: maximal number of bases next to position of interest that are overlaped by a single read
#' - `a`: reads mapping to sequence left of the position of interest
#' - `b`: reads mapping to sequence right of the position of interest
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


#' Map the re-quantification result on the junction-transcript centric data by hash id.
#'
#' @param path_to_easyquant_folder The path to easyquant folder
#' @param junc_tib The junction-transcript centric data predicted from WES data. A unique context sequence needs to be described by a hash id in a column `hash_id`
#'
#' @return Extended junction tibble with re-quantification results.  The ollowing columns are added:
#' -  `pos`: position of interest relative to input sequence
#' - `junc`: reads overlapping the position of interest
#' - `span`: read pairs spanning the position of interest
#' - `anch`: maximal number of bases next to position of interest that are overlaped by a single read
#' - `a`: reads mapping to sequence left of the position of interest
#' - `b`: reads mapping to sequence right of the position of interest
#
#'
#' @examples
#'
#'
#'@import dplyr
#'
#'@export
map_requant <- function(path_to_easyquant_folder, junc_tib) {
  dat_easyqant <- read_requant(path_to_easyquant_folder)
  dat_junc <- junc_tib %>%
    left_join(dat_easyqant, by = c("cts_id" = "name"))
  return(dat_junc)
}


