#' Imports re-quantification results with easyquant
#'
#' @param path The path to easyquant folder
#'
#' @return A tibble with with the requantification results. One junction might be re-presented by multiple context sequences, i.e. transcripts
#'
#' @examples
#'
#'
#'@import dplyr
#'
#'@export
read_requant <- function(path_folder){
  path_folder <- "/scratch/info/projects/SUMMIT/WP1.1/alternative_splicing/data/polyA/20210217_easyquant/S01_rep1"
  path.to.easyquant.file <- paste(path_folder, "quantification.tsv" ,sep = "/" )
  if(!file.exists(path.to.easyquant.file)){
    stop("quantification.tsv file is missing")
  }
  dat.easyqant <- path.to.easyquant.file %>%
    readr::read_delim(delim = "\t") %>%
    tidyr::separate(name, into = c("chrom", "start", "end", "strand", "transcript_id"), sep = "_", remove = F) %>%
    mutate(junc_id = paste(chrom, start, end, strand ,sep = "_")) %>%
  return(dat.easyqant)
}


#' Transforms re-quantification results into junc_id centric tibble
#'
#' @param dat The path to easyquant folder
#'
#' @return A tibble with with the requantification results. One junction is represented by one row.
#'
#' @examples
#'
#'
#'@import dplyr
#'
#'@export
transform_requant <- function(dat) {
  dat_junc <- dat %>%
    group_by(junc_id) %>%
    summarise(
      mean_junc =  mean(junc),
      max_junc = max(junc),
      mean_span =  mean(span),
      max_span = max(span),
      number_transcripts = length(junc),
      transcripts = paste(transcript_id, collapse = ",")
    )


}
