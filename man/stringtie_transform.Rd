% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_stringtie.R
\name{stringtie_transform}
\alias{stringtie_transform}
\title{Imports StringTie assembled transcripts and transforms the raw output
into standardized junction output format}
\usage{
stringtie_transform(gtf.file)
}
\arguments{
\item{gtf.file}{The path to StringTie GTF file}
}
\value{
A tibble in standardized junction format
}
\description{
\itemize{
\item GitHub: https://github.com/gpertea/stringtie
\item Paper: https://doi.org/10.1038/nbt.3122
}
}
