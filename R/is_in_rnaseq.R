#' test if junction was found in corresponding RNA-seq data
#'
#' @param junc_id vector of junction id to test
#' @param rna_juncs a vector of junctions that were identified in RNA-seq data
#'
#' @return logical vector of same length as `junc_id` indicating
#' if the input `junc_id` was found in RNA-seq data of the same sample
#'
#' @export
is_in_rnaseq <- function(junc_id, rna_juncs){

  in_rna_junc <- junc_id %in% rna_juncs

  return(in_rna_junc)

}
