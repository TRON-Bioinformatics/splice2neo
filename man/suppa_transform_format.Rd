% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_suppa2.R
\name{suppa_transform_format}
\alias{suppa_transform_format}
\title{Transforms SUPPA2 ioe event files file into standardized junction format.}
\usage{
suppa_transform_format(tib_list)
}
\arguments{
\item{tib_list}{A list with tibbles that contain the SUPPA2 output - each for one type of alternative splicing.
These types can be "A5SS", "A3SS", "cassette_exon", "intron_retention", "mutex_exons".}
}
\value{
A tibble in standardized junction format, combining all alternative
splicing classes that are were determined with SUPPA2
}
\description{
Transforms SUPPA2 ioe event files file into standardized junction format.
}
