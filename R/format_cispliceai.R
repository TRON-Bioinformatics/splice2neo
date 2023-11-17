
#' Formats CI-SpliceAI with -t flag output and filter for predicted effects
#'
#' Reformat the data for each annotated effect per row, filters effects
#' to have a probability score not NA and score > 0, and removes gene symbol from
#' data to make non-redundant output.
#'
#' @param spliceai_variants [tibble][tibble::tibble-package] with parsed
#' spliceAI mutations from \code{\link{parse_spliceai}}
#'
#' @return A [tibble][tibble::tibble-package] with splicing effects per row
#'
#' @examples
#' cispliceai_file <- system.file("extdata", "cispliceai_thresh_output.vcf", package = "splice2neo")
#' df <- parse_cispliceai_thresh(cispliceai_file)
#' format_cispliceai_thresh(df)
#'
#' @seealso \code{\link{parse_cispliceai_thresh}}, \code{\link{annotate_mut_effect}}
#' @export
format_cispliceai_thresh <- function(cispliceai_variants){
  
  #format columns
  cispliceai_variants <- cispliceai_variants %>%
    mutate(score = as.numeric(score),
           pos_rel = as.integer(pos_rel), 
           effect = as.factor(effect)) %>%
    
    # filter effects without probability score given
    filter(!is.na(score) & score > 0) %>%
    
    # add unique IDs for mutation
    mutate(
      mut_id = str_c(CHROM, POS, REF, ALT, sep = "_"),
      chr = CHROM,
      pos = as.integer(POS) + pos_rel
    ) %>%
    
    # keep only relevant columns
    dplyr::select(mut_id, effect, score, chr, pos_rel, pos) %>%
    dplyr::distinct()
}