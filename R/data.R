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


#' A tibble in junction format but unsorted columns
#'
#' @format A tibble frame with 1 and 8 variables:
#'
#' @source random
"unsorted_junc_df"

#' A tibble containing events from alternative 3' splice sites identified by
#' Spladder
#'
#' @format A tibble frame with 1 and 15 variables:
#'
#' @source random
"spladder_output.a3ss"

#' A tibble containing events from alternative 5' splice sites identified by
#' Spladder
#'
#' @format A tibble frame with 1 and 15 variables:
#'
#' @source random
"spladder_output.a5ss"

#' A tibble containing events from exon skipping identified by
#' Spladder
#'
#' @format A tibble frame with 1 and 15 variables:
#'
#' @source random
"spladder_output.exonskip"


#' A tibble containing events from intron retention identified by
#' Spladder
#'
#' @format A tibble frame with 1 and 15 variables:
#'
#' @source random
"spladder_output.intronreten"

#' A tibble containing events from mutually exclusive exons identified by
#' Spladder
#'
#' @format A tibble frame with 1 row and 17 variables:
#'
#' @source random
"spladder_output.mutexon"

#' A list of tibbles with each list containing events from a different
#' splicing class identified by Spladder
#'
#' @format A list frame with tibbles called: A3SS, A5SS, casette_exon,
#' intron_retention and mutex_exons
#'
#' @source random
"spladder_output"


#' A tibble with junction in standardized format. example data
#'
#' @format A tibble frame with 11 variables
#'
#' @source random
"leafcutter_juncs"

#' A tibble with junction in standardized format.  example data
#'
#' @format A tibble frame with 11 variables
#'
#' @source random
"spladder_juncs"


#' A tibble with canonical junctions and their source (comma, separated).
#' example data
#'
#' @format A tibble frame with 2 variables
#'
#' @source random
"canonical_juncs"

