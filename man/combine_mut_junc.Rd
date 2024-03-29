% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine_mut_junc.R
\name{combine_mut_junc}
\alias{combine_mut_junc}
\title{Combines data sets with junctions from several different sources}
\usage{
combine_mut_junc(junc_data_list)
}
\arguments{
\item{junc_data_list}{A named list of tibbles with junctions. the name
should be the source, e.g. the tool name, such as \code{spliceai} or \code{pangolin}.
The individual tibbles should at least contain the columns \code{mut_id},
\code{tx_id}, \code{junc_id} and might have individual sets of other tool/source
specific columns. It is recommended that these input tables are unique
with respect to \code{mut_id}, \code{tx_id}, \code{junc_id}, e.g. by using \code{unique_mut_junc()} or
\code{unique_junc_mmsplice()}.}
}
\value{
A combined data set with unique junctions based on the columns
\code{mut_id}, \code{tx_id}, \code{junc_id}.
Additional columns in the input data.frames will be prefixed with the
tool/source name followed by an underscore \verb{_}.
E.g. the column \code{score} in the input data sets \code{spliceai} becomes
\code{spliceai_score}. Additionally, for each tool/source with name \verb{<name>}
an additional column \verb{<name>_detected} will be added to indicate if a
given junction was detected with the indicated tool.
}
\description{
Combines data sets with junctions from several different sources
}
