
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

  gr <- rtracklayer::import.bed(bed_file)
  chr <- GenomeInfoDb::seqnames(gr)
  strand <- as.character(BiocGenerics::strand(gr))

  if(type == "exon-exon"){

    start <-BiocGenerics::start(gr)
    end <- BiocGenerics::end(gr)

  } else if(type == "intron"){
    start <-BiocGenerics::start(gr) - 1
    end <- BiocGenerics::end(gr) + 1
  }

  junc_id <- generate_junction_id(chr,start, end, strand)

  return(junc_id)
}
