#' An example transcript annotation GRangesList from a subset of human genome
#'
#' This is a GRangesList object with transcripts consisting of exon ranges.
#' This data set is a subset of annotations of human (hg19) restricted
#' to this regions:
#'
#'  - chr2:152000000-180000000
#'  - chr17:41100000-41280000
#'
#' @format A GRangesList object with GRanges objects.
#' @source \url{ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz}
"toy_transcripts"

#' An example dataset of full transcript ranges from a subset of human genome.
#'
#' This is a GRanges object with full transcript ranges.
#' This data set is a subset of annotations of human genome (hg19)
#' restricted to this regions:
#'
#'  - chr2:152000000-180000000
#'  - chr17:41100000-41280000
#'
#' @format A GRanges object with regions.
#' @source \url{ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz}
"toy_transcripts_gr"
