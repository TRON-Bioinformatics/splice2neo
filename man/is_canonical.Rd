% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_canonical.R
\name{is_canonical}
\alias{is_canonical}
\title{test if junction is canonical junction}
\usage{
is_canonical(junc_id, ref_junc, exons_gr)
}
\arguments{
\item{junc_id}{vector of junction id to test}

\item{ref_junc}{vector of canonical reference junctions}

\item{exons_gr}{GRanges of canonical exons}
}
\value{
logical vector of same length as \code{junc_id} indicating
if the input \code{junc_id} is in the set of reference junctions \code{ref_junc} or
consists of two directly adjacent positions (intron retention) and
overlaps (completely within) a canonical exon from \code{exons_gr}.
}
\description{
test if junction is canonical junction
}
