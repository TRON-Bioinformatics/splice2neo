# STAR functions ----------------------------------------------------------

#' Imports tabular STAR SJ.out.tab splice junctions
#'
#' @param file.name Path to STAR SJ.out.tab
#'
#' @return A tibble with one splice junction per row and columns
#'  chromosome, intron_start, intron_end, strand_enc, motif_enc,
#'  annotation, uniquely_mapping_reads, multi_mapping_reads, alignment_overhang
#'
#' @examples
#'
#' path <-  system.file("extdata", "test_star_SJ.out.tab", package = "splice2neo")
#' splice2neo:::import_star_sj(path)
#'
#' @import readr
#' @keywords internal
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
  junc.file <- read_tsv(
    file.name,
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
#' Junctions without strand annotation "*" are filtered out.
#'
#' @examples
#' path <-  system.file("extdata", "test_star_SJ.out.tab", package = "splice2neo")
#' star_raw <- splice2neo:::import_star_sj(path)
#' splice2neo:::transform_star_sj(star_raw)
#'
#' @import dplyr
#' @keywords internal
transform_star_sj <- function(tib){

  tib <- tib %>%
    mutate(
      junction_start = as.character(intron_start - 1),
      junction_end = as.character(intron_end + 1),
      strand = case_when(
        strand_enc == "0"  ~ "*",
        strand_enc == "1"  ~ "+",
        strand_enc == "2"  ~ "-",
      ),
      Gene = NA,
      class = NA
    ) %>%
    filter(strand != "*") %>%
    mutate(
      junc_id = generate_junction_id(chromosome, junction_start, junction_end, strand)
    ) %>%
    select(
      junction_start,
      junction_end,
      strand,
      chromosome,
      Gene,
      class,
      junc_id,
      uniquely_mapping_reads,
      multi_mapping_reads
    )

  return(tib)
}

#' Imports "*SJ.out.tab" from STAR and transforms the raw output
#' into standardized junction output format.
#'
#' @param path The path to STAR SJ.out.tab
#'
#' @return A tibble in standardized junction format
#'
#' @examples
#'
#' path <-  system.file("extdata", "test_star_SJ.out.tab", package = "splice2neo")
#' parse_star_sj(path)
#'
#' @import readr
#' @export
parse_star_sj <- function(path) {

  dat.junc <- path %>%
    import_star_sj() %>%
    transform_star_sj()
  return(dat.junc)
}
