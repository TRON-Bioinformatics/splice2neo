% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_star.R
\name{transform_star_sj}
\alias{transform_star_sj}
\title{Transforms STAR intermediate table into standardized junction format}
\usage{
transform_star_sj(tib, keep_unstranded = FALSE)
}
\arguments{
\item{tib}{STAR *SJ.out.tab as tibble}

\item{keep_unstranded}{Logical, to keep junctions without strand annotation.
If FALSE (default) junctions without strand annotation "*" are filtered out,
otherwise they are considered as two separate junctions with "+" and "-" strand.}
}
\value{
A tibble in standardized junction format
}
\description{
Transforms STAR intermediate table into standardized junction format
}
\examples{
path <-  system.file("extdata", "test_star_SJ.out.tab", package = "splice2neo")
star_raw <- splice2neo:::import_star_sj(path)
splice2neo:::transform_star_sj(star_raw)

}
\keyword{internal}
