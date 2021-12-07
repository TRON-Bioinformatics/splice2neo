
#' Transform a bed file into junction format
#'
#' @param bed_file Path to the bed file
#' @param type Use for this parameter
#'   -  `exon-exon` if the bed file defines exon-exon boundaries
#'   -  `intron` if the bed file defines introns
#'
#' @return A character vector with junction IDs
#'
#'
#'
#' @import dplyr
#' @export
bed_to_junc <- function(bed_file, type = "exon-exon"){

  stopifnot(is.character(bed_file))
  stopifnot(type %in% c("exon-exon", "intron"))

  if(type == "exon-exon"){
    gr <- rtracklayer::import.bed(bed_file)
    junc_id <- paste(GenomeInfoDb::seqnames(gr), BiocGenerics::start(gr), BiocGenerics::end(gr), BiocGenerics::strand(gr), sep = "_")

  } else if(type == "intron"){
    gr <- rtracklayer::import.bed(bed_file)
    junc_id <- paste(GenomeInfoDb::seqnames(gr), BiocGenerics::start(gr) - 1, BiocGenerics::end(gr) + 1, BiocGenerics::strand(gr), sep = "_")

  }

  return(junc_id)
}
