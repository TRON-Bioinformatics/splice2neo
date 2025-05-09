% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annotate_jet.R
\name{annotate_potential_jet}
\alias{annotate_potential_jet}
\title{Annotate if splice sites of a junction overlap with genomic positions of transposable elements.}
\usage{
annotate_potential_jet(df, rmsk)
}
\arguments{
\item{df}{A data.frame with splice junctions in rows and at least the column:
\itemize{
\item \code{junc_id} junction id consisting of genomic coordinates
}}

\item{rmsk}{\code{\link[GenomicRanges]{GRanges}} of transposable elements (e.g. RepeatMasker)}
}
\value{
A data.frame as the input with additional columns annotating overlaps with transposable elements.
\itemize{
\item \code{potential_jet} Indicator if junction is potentially a JET
\item \code{left_side_retroelement}: Name of transposable element overlapping left splice site
\item \code{left_side_retroelement_class}: Class of transposable element overlapping left splice site
\item \code{right_side_retroelement}: Name of transposable element overlapping right splice site
\item \code{right_side_retroelement_class}: Class of transposable element overlapping right splice site
}
}
\description{
Annotate if splice sites of a junction overlap with genomic positions of transposable elements.
}
\examples{


rmsk <- readr::read_tsv(system.file("extdata", "rmsk_hg19_subset.tsv.gz", package = "splice2neo"))
rmsk <- GenomicRanges::makeGRangesFromDataFrame(rmsk, keep.extra.columns = TRUE)
junc_df <- tibble::tibble(
 junc_id = c("chr2:152407458-152408252:-")
)

annotate_potential_jet(junc_df, rmsk)

}
