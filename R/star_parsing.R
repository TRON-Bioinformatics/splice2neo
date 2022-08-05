# STAR functions ----------------------------------------------------------

#' Imports tabular STAR SJ.out.tab splice junctions
#'
#' @param file.name Path to STAR SJ.out.tab
#'
#' @return A tibble with one splice junction per row and columns
#'  chromosome, intron_start, intron_end, strand_enc, motif_enc, 
#'  annotation, uniquely_mapping_reads, multi_mapping_reads, alignment_overhang
#'
#' @import readr
#' @export
import_star_sj <- function(file.name) {
  if(!file.exists(file.name)){
    stop("STAR SJ.out.tab is missing")
  }
  star.cols <-
    c('chromosome',
      'intron_start',
      'intron_end',
      'strand_enc',
      'motif_enc',
      'annotation_enc',
      'uniquely_mapping_reads',
      'multi_mapping_reads',
      'alignment_overhang'
    )
  junc.file <- read_delim(
    file.name,
    delim = '\t',
    col_names = star.cols,
    col_types = readr::cols(
      .default = readr::col_integer(),
      chromosome = readr::col_character()
    )
  )

  return(junc.file)
}


#' Transforms STAR intermediate table into standardized junction format
#'
#' @param tib STAR *SJ.out.tab as tibble
#'
#' @return A tibble in standardized junction format
#'
#'
#' @import dplyr
#' @export
transform_star_sj <- function(tib){
  strand_info <-
    c(
      "1" = "+",
      "2" = "-",
      "0" = "*"
    )
  tib <- tib %>%
    mutate(
      junction_start = as.character(intron_start - 1),
      junction_end = as.character(intron_end + 1),
      strand = strand_info[as.character(strand_enc)],
      Gene = NA,
      class = NA,
      AS_event_ID = NA
    ) %>%
    filter(strand != "*") %>%
    mutate(
      junc_id = generate_junction_id(chromosome, junction_start, junction_end, strand)
    ) %>%
    sort_columns()
  return(tib)
}

#' Imports "*SJ.out.tab" from STAR and transforms the raw output
#' into standardized junction output format.
#'
#' @param path The path to STAR SJ.out.tab
#'
#' @return A tibble in standardized junction format
#' @examples
#' path <-  system.file("extdata", "test_star_SJ.out.tab", package = "splice2neo")
#' star_juncs <- star_transform(path)
#' star_juncs
#' 
#' @import readr
#' @export
star_transform <- function(path) {
  file.junc <- path
  dat.junc <- file.junc %>%
    import_star_sj() %>%
    transform_star_sj()
  return(dat.junc)
}
