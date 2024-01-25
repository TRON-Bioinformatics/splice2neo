# IRFinder functions ----------------------------------------------------------

#' Imports tabular IRFinder retained introns predictions
#'
#' @param file.name Path to "IR-(non)dir-(val).txt"
#'
#' @return A tibble with one retained intron call per row
#'
#' @import readr
import_irfinder_txt <- function(file.name) {
  if(!file.exists(file.name)){
    stop("IRFinder-IR-nondir.txt file is missing...")
  }
  intron.file <- read_tsv(
    file.name,
    col_types = readr::cols(
      Chr = readr::col_character(),
      Start = readr::col_integer(),
      End = readr::col_integer()
    )
  )

  return(intron.file)
}


#' Transforms IRFinder intermediate table into standardized junction format
#'
#' @param tib IRFinder IR-(non)dir-(val).txt as tibble
#'
#' @return A tibble in standardized junction format
#'
#'
#' @import dplyr
transform_irfinder_txt <- function(tib){
  tib <- tib %>%
    mutate(
      junc_id1 = generate_junction_id(Chr, Start, Start + 1, Strand),
      junc_id2 = generate_junction_id(Chr, End, End + 1, Strand),
      junc_id = paste(junc_id1, junc_id2, sep = ";"),
      strand = Strand,
      Gene = stringr::str_split(Name, "/") %>% purrr::map_chr(., 2),
      overlapping_feature = stringr::str_split(Name, "/") %>% purrr::map_chr(., 3),
      class = 'intron_retention',
      AS_event_ID = Name
    ) %>%
    tidyr::separate_rows(junc_id, sep=";") %>%
    tidyr::separate(
      junc_id,
      into = c("chromosome", "junction_start_end", "strand"),
      sep = ":",
      remove = F
    ) %>%
    tidyr::separate(
      junction_start_end,
      into = c("junction_start", "junction_end"),
      sep = "-",
      remove = T
    ) %>%
    dplyr::select(
      chromosome, 
      junction_start,
      junction_end,
      strand,
      junc_id,
      Gene,
      class,
      AS_event_ID,
      overlapping_feature,
      IRratio,
      Warnings,
      IntronDepth
    )
  return(tib)
}

#' Filter IRFinder intermediate table to remove likely 
#' false postive IR predictions  
#'
#' @param tib IRFinder intermediate tibble with columns Warnings and IRratio
#' @param warnings Filter out IR calls with Warnings
#' @param ratio_cutoff Filter out IR calls with IRRatio <= this parameter
#'
#' @return A tibble in standardized junction format
#'
#'
#' @import dplyr
#' @export
filter_irfinder_txt <- function(tib, warnings=FALSE, ratio_cutoff=0.1){
  stopifnot("Argument ratio_cutoff must be a value between 0 and 1!" = ratio_cutoff > 0 & ratio_cutoff <= 1)
  filtered_introns <- tib %>%
    dplyr::filter(IRratio >= ratio_cutoff)
  
  if (isFALSE(warnings)){
    message('Removing IR calls with warnings...')
    filtered_introns <- filtered_introns %>%
      dplyr::filter(Warnings == '-')
  }
  
  message('Removing IR calls with overlapping features...')
  filtered_introns <- filtered_introns %>%
    dplyr::filter((stringr::str_split(Name, "/") %>% purrr::map_chr(., 3)) == 'clean')

  return(filtered_introns)

}

#' Imports "IRFinder-IR-nondir.txt" from IRFinder short mode
#' and transforms the raw output into standardized junction output format.
#'
#' @param path The path to "IR-(non)dir-(val).txt"
#' @param warnings Filter out IR calls with Warnings
#' @param irratio  Filter out IR calls with IRRatio <= this parameter
#' @param cnn Use CNN validated table as input instead of standard output table
#'
#' @return A tibble in standardized junction format
#' 
#' path <- system.file("extdata", "", package = "splice2neo")
#' ir_juncs <- parse_irfinder_txt(path)
#' ir_juncs
#' 
#' @import readr fs
#' @export
parse_irfinder_txt <- function(path, warnings=FALSE, irratio=0.1, cnn=FALSE) {
  if (isTRUE(cnn)) {
    message('Importing CNN validated IR predictions...')
    file.intron <- fs::path_join(c(path, 'IRFinder-IR-nondir-val.txt'))  
  } else {
    message('Importing standard non-directional IR predictions...')
    file.intron <- fs::path_join(c(path, 'IRFinder-IR-nondir.txt'))
  }
  dat.junc <- file.intron %>%
    import_irfinder_txt() %>%
    filter_irfinder_txt(warnings = warnings, ratio_cutoff = irratio) %>%
    transform_irfinder_txt()
  return(dat.junc)
}
