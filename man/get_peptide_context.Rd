% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_peptide.R
\name{get_peptide_context}
\alias{get_peptide_context}
\title{Get peptide context sequence given a tibble with
\code{intron_retention}
\code{normalized_protein_junc_pos}
\code{is_first_reading_frame}
\code{frame_shift}
\code{exon1_end_AA}
\code{exon1_end_AA_WT}
\code{exon2_start_AA_WT},
\code{protein_len}
\code{protein_len}
\code{protein_length_difference}
\code{protein}}
\usage{
get_peptide_context(df, flanking_size = 14)
}
\arguments{
\item{df}{Data frame with information about position of junction etc.}

\item{flanking_size}{number amino acids left and right of the breakpoint or novel sequence part}
}
\description{
Get peptide context sequence given a tibble with
\code{intron_retention}
\code{normalized_protein_junc_pos}
\code{is_first_reading_frame}
\code{frame_shift}
\code{exon1_end_AA}
\code{exon1_end_AA_WT}
\code{exon2_start_AA_WT},
\code{protein_len}
\code{protein_len}
\code{protein_length_difference}
\code{protein}
}
\keyword{internal}
