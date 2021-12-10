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
#' - `name`: name of the input sequence. this will be a hash id.
#' -  `pos`: position of interest relative to input sequence
#' - `junc`: reads overlapping the position of interest
#' - `span`: read pairs spanning the position of interest
#' - `anch`: maximal number of bases next to position of interest that are overlapped by a single read
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


#' Maps the re-quantification result from EasyQuant. on the junction-transcript centric tibble by hash id.
#'
#' @param path_to_easyquant_folder The path to EasyQuant folder
#' @param junc_tib The junction-transcript centric tibble, i.e. each row represents an altered transcript. Must contain a column `hash_id` with hash ids that relate to the column `name` in the Easyquant output.
#'
#' @return Extended junction-transcript tibble with re-quantification results.  The following columns are added:
#' -  `pos`: position of interest relative to input sequence
#' - `junc`: reads overlapping the position of interest
#' - `span`: read pairs spanning the position of interest
#' - `anch`: maximal number of bases next to position of interest that are overlapped by a single read
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


