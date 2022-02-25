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

#' An example CDS annotation GRangesList from a subset of human genome
#'
#' This is a GRangesList object with GenomicRanges for each transcript consisting of CDS ranges.
#' This data set is a subset of annotations of human (hg19) restricted
#' to this regions:
#'
#'  - chr2:152000000-180000000
#'  - chr17:41100000-41280000
#'
#' @format A GRangesList object with GRanges objects.
#' @source \url{ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz}
"toy_cds"

#' An example dataset of 18 splice junctions in `junc_id` format from human
#' (hg19) genome as character vector
#'
#' @format A character vector with `junc_id`s in the format:
#'  `<chr>_<pos1>_<pos2>_<strand>`
#'
#' @source They were created from the spliceAI example file
"toy_junc_id"

#' An example dataset of 18 transcript IDs matching to the junctions in `toy_junc_id`
#'
#' @format A character vector with ENSEMBL transcript ids
#'
#' @source They were created using transcript annotations from
#' \url{https://www.ensembl.org}
#'
"toy_junc_id_enst"

#' An example dataset of 18 splice junctions from human
#' (hg19) genome as a data.frame
#'
#' The data.frame has the following columns:
#'
#'   - `junc_id` see \link{toy_junc_id}
#'   - `tx_id` see \link{toy_junc_id_enst}
#'
#' @source They were created from the spliceAI example file
"toy_junc_df"

#' A tibble in junction format but unsorted columns
#'
#' @format A tibble frame with 1 and 8 variables:
#'
#' @source random
"unsorted_junc_df"


#' A tibble with canonical junctions and their source (comma, separated).
#' example data
#'
#' @format A tibble frame with 2 variables
#'
#' @source random
"canonical_juncs"


