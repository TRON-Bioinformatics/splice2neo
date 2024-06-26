% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_gtf.R
\name{parse_gtf}
\alias{parse_gtf}
\title{Parse a GFF/GTF file as \code{\link[GenomicRanges]{GRangesList}} of exons}
\usage{
parse_gtf(file, ...)
}
\arguments{
\item{file}{path to a GFF or GTF file (See \code{\link[GenomicFeatures]{makeTxDbFromGFF}}))}

\item{...}{other options passed to file path to \code{\link[GenomicFeatures]{makeTxDbFromGFF}})}
}
\value{
A \code{\link[GenomicRanges]{GRangesList}} of exons grouped by transcripts
}
\description{
Parse a GFF/GTF file as \code{\link[GenomicRanges]{GRangesList}} of exons
}
\examples{
gff_file <- system.file("extdata","GFF3_files","a.gff3",package="splice2neo")

parse_gtf(gff_file)

}
