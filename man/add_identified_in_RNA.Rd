% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardize_junctions_RNA_tools.R
\name{add_identified_in_RNA}
\alias{add_identified_in_RNA}
\title{This is a wrapper function to directly map the information if a junction predicted from WES data was found in RNA-seq by Regtools or SplAdder}
\usage{
add_identified_in_RNA(mutation_juncs, path_to_spladder, path_to_regtools)
}
\arguments{
\item{mutation_juncs}{The junction-transcript centric data predicted from WES data.}

\item{path_to_spladder}{The path to the results from RNA-seq analysis with SplAdder}

\item{path_to_regtools}{The path to the results from RNA-seq analysis with Regtools}
}
\value{
The junction-transcript centric data predicted from WES data is extended by the information if a respective aberrant junctions was
identified by spladder or regtools (\code{identified_by_regtools}, \code{identified_by_spladder})
}
\description{
This is a wrapper function to directly map the information if a junction predicted from WES data was found in RNA-seq by Regtools or SplAdder
}
