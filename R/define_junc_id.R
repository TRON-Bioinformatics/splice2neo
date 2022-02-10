#' Given the chromosome, junction start, junction end and strand,
#' a junction id is created that follows the format: `<chr>:<start>-<end>:<strand>`
#'
#' @param chr The chromosome
#' @param start The start of the junction
#' @param end The end of the junction
#' @return The junction id in the format `<chr>:<start>-<end>:<strand>`
#'
#' @examples
#'
#' generate_junction_id(chr = "chr1", start = 50, end = 100, strand = "+")
#'
#'@import dplyr
#'
generate_junction_id <- function(chr, start, end, strand){

  stopifnot(strand %in% c("+", "-"))

  # create junction id
  junc_id <- str_c(chr, ":", as.integer(start),"-", as.integer(end),":", strand)

  return(junc_id)

}

#' Given the chromosome, junction start, junction end and strand,
#' a junction id is created that follows the format: `<chr>:<start>-<end>:<strand>`
#'
#' @param junc_id The junction id in the format `<chr>:<start>-<end>:<strand>`
#' @return The breakpoint kit id in the format `<chr>:<start>-<chr>:<end>:<strand>`
#'
#' @examples
#'
#'
#'
#'@import dplyr
#'
junc2breakpoint <- function(junc_id){

}
