#' Given the chromosome, junction start, junction end and strand,
#' a junction id is created that follows the format: `<chr>:<start>-<end>:<strand>`
#'
#' @param chr The chromosome
#' @param start The start of the junction
#' @param end The end of the junction
#' @param strand Strand information. "+" or "-"
#' @return The junction id in the format `<chr>:<start>-<end>:<strand>`
#'
#' @examples
#'generate_junction_id(chr = "chr1", start = 50, end = 100, strand = "+")
#'
#'@import dplyr
#'@export
generate_junction_id <- function(chr, start, end, strand){

  stopifnot(strand %in% c("+", "-"))

  # create junction id
  junc_id <- str_c(chr, ":", as.integer(start),"-", as.integer(end),":", strand)

  return(junc_id)

}

#' Transforms the junction id into the breakpoint id of the format `<chr>:<start>-<chr>:<end>`
#'
#' @param junc_id The junction id in the format `<chr>:<start>-<end>:<strand>`
#' @return The breakpoint kit id in the format `<chr>:<start>-<chr>:<end>`
#'
#' @examples
#' junc_id <- generate_junction_id(chr = "chr1", start = 50, end = 100, strand = "+")
#' junc2breakpoint(junc_id)
#'
#'@import dplyr
#'@export
junc2breakpoint <- function(junc_id){

  # split junction id in chr, start, end
  junc_split <- str_split_fixed(junc_id, ":", n = 3)
  chr <- junc_split[,1]
  start_end <- str_split_fixed(junc_split[,2], "-", n = 2)
  start <- start_end[,1]
  end <- start_end[,2]

  breakpoint_id <- stringr::str_c(chr,":", start, "-", chr, ":", end)


}

#' Given the strand Transforms the breakpoint id into the junction id of the format `<chr>:<start>-:<end>:<strand>`
#'
#' @param breakpoint_id The junction id in the format `<chr>:<start>-<chr>:<end>`
#' @param strand Strand information. "+" or "-"
#' @return The junction id in the format `<chr>:<start>-<end>:<strand>`
#'
#' @examples
#'
#' breakpoint_id <- "chr1:500000-chr1:1000000"
#' breakpoint2junc(breakpoint_id = breakpoint_id, strand = "+")
#'
#'@import dplyr
#'@export
breakpoint2junc <- function(breakpoint_id, strand){

  stopifnot(strand %in% c("+", "-"))

  # split junction id in chr, start, end
  id_split <- str_split_fixed(breakpoint_id, "-", n = 2)
  chr_1 <- str_split_fixed(id_split[,1], ":", n = 2)[,1]
  start <- str_split_fixed(id_split[,1], ":", n = 2)[,2]
  chr_2 <- str_split_fixed(id_split[,2], ":", n = 2)[,1]
  end <- str_split_fixed(id_split[,2], ":", n = 2)[,2]

  stopifnot("the breakpoint id can only be converted into a junction id, if start and end are located on the same chromosome"=chr_1 == chr_2)

  junc_id <- generate_junction_id(chr_1, start, end, strand)


}


