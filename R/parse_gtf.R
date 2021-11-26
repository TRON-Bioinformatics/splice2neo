

#' Parse a GFF/GTF file as \code{\link[GenomicRanges]{GRangesList}} of exons
#'
#' @param file path to a GFF or GTF file (See \code{\link[GenomicFeatures]{makeTxDbFromGFF}}))
#' @return A \code{\link[GenomicRanges]{GRangesList}} of exons grouped by transcripts
#'
#' @examples
#' gff_file <- system.file("extdata","GFF3_files","a.gff3",package="GenomicFeatures")
#'
#' parse_gtf(gff_file)
#'
#' @export
parse_gtf <- function(file){

  txdb <- suppressWarnings(GenomicFeatures::makeTxDbFromGFF(file))
  transcripts <- GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)

  return(transcripts)
}
