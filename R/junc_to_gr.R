

#' Convert a splice junction ID into a GRanges object
#'
#' @param junc_id A character vector with a junction IDs
#'
#' @return A \code{\link[GenomicRanges]{GRanges}} object
#'
#' @examples
#' junc_to_gr("chr1_5_10_+")
#'
#' junc_id <- c("chr1_5_10_+", "chr_special1_5_20_-")
#' junc_to_gr(junc_id)
#'
#' @export
junc_to_gr <- function(junc_id){

  stopifnot(!is.na(junc_id))

  fields <- stringr::str_match(junc_id, "(.*)_(\\d+)_(\\d+)_([*+-])")

  chr <- fields[,2]
  start <- fields[,3] %>% as.integer()
  end <- fields[,4] %>% as.integer
  strand <- fields[,5]

  gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(start, end), strand = strand,
                         names = junc_id)
  names(gr) <- junc_id
  return(gr)
}
