

#' Convert a splice junction ID into a GRanges object
#'
#' @param junc_id A character vector with a junction IDs
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object
#'
#' @examples
#' junc_to_gr("chr5:10-50:+")
#'
#' junc_id <- c("chr5:10-50:+", "chr_special1:5-20:-")
#' junc_to_gr(junc_id)
#'
#' @export
junc_to_gr <- function(junc_id){

  stopifnot(!is.na(junc_id))

  gr <- GenomicRanges::GRanges(junc_id)

  return(gr)
}
