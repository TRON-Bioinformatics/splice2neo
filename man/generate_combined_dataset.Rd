% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/standardize_junctions_RNA_tools.R
\name{generate_combined_dataset}
\alias{generate_combined_dataset}
\title{Combines tibbles with junctions from any number of RNA-seq tools into a
combined dataset of expressed splice junctions}
\usage{
generate_combined_dataset(rna_junc_data_list)
}
\arguments{
\item{rna_junc_data_list}{A named list with junction tibbles in standardized
format.}
}
\value{
A combined table with unique junctions. The columns
identified_by_\{name\} contains information which tools identified the given
junction
}
\description{
Combines tibbles with junctions from any number of RNA-seq tools into a
combined dataset of expressed splice junctions
}
\examples{
path <-  system.file("extdata", "", package = "splice2neo")
spladder_juncs <- spladder_transform(path)
path <-  system.file("extdata", "test_regtools_Aligned.out.sorted.bam.junc", package = "splice2neo")
regtools_juncs <- regtools_transform(path)
dat.combined <- generate_combined_dataset(list("spladder" = spladder_juncs,
                                               "regtools" = regtools_juncs))

}
