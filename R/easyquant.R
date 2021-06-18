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
  path_folder <- "/scratch/info/projects/SUMMIT/WP1.1/alternative_splicing/data/polyA/20210217_easyquant/S01_rep1"
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
#' @param dat_mutation The junction-transcript centric data predicted from WES data. A unique context sequence needs to be described by a hash id in a column `hash_id`
#'
#' @return Extended junction tibble with re-quantification results.  Followin columns are added:
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
map_requant <- function(path_to_easyquant_folder, dat_mutation) {
  dat_easyqant <- read_requant(path_to_easyquant_folder)
  dat_junc <- dat_mutation %>%
    left_join(dat_easyqant, by = c("hash_id" = "name"))
  return(dat_junc)
}


