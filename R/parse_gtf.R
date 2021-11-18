

#' Parse a GFF/GTF file as \code{\link[GenomicRanges]{GRangesList}} of exons
#'
#' @param file path to a GFF or GTF file (See \code{\link[makeTxDbFromGFF]{GenomicFeatures}}))
#' @return A \code{\link[GenomicRanges]{GRangesList}} of exons grouped by transcripts
#'
#' @export
parse_gtf <- function(file){

  txdb <- GenomicFeatures::makeTxDbFromGFF(file)
  transcripts <- GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)

  return(transcripts)

}

canonical_junctions <- function(grl){

}
