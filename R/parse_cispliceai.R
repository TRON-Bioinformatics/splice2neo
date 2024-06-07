
#' Parse VCF output file from CI-SpliceAI with -t flag as table
#'
#' requires the vcfR package.
#'
#' @param vcf_file path to a single VCF file, the output of spliceAI
#' @return a [tibble][tibble::tibble-package] with one row per variant
#' annotation. Each input variant can have multiple annotations.
#'
#' @examples
#'
#' cispliceai_file <- system.file("extdata", "cispliceai_thresh_output.vcf", package = "splice2neo")
#' parse_cispliceai_thresh(cispliceai_file)
#'
#' @seealso \code{\link{format_cispliceai_thresh}}, \code{\link{annotate_mut_effect}}
#' @import dplyr stringr tidyr
#' @export
parse_cispliceai_thresh <- function(vcf_file){

  # read in vcf file
  vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)

  # get fixed fields for variants
  if (nrow(vcf) > 1){
    fix_df <-  vcfR::getFIX(vcf) %>%
      tibble::as_tibble() %>%
      mutate(Key = row_number())
  } else {
    fix_df <-  vcfR::getFIX(vcf) %>%
      tibble::as_tibble_row() %>%
      mutate(Key = row_number())
  }

  # get annotations
  splice_df <- vcf %>%
    vcfR::extract_info_tidy(info_fields = c("CI-SpliceAI")) %>%
    mutate(fields = `CI-SpliceAI` %>% stringr::str_split(",")) %>%
    dplyr::select(-`CI-SpliceAI`) %>%
    tidyr::unnest(fields) %>%
    dplyr::filter(!is.na(fields)) %>%

  # parse multiple splice effects from fields column
    mutate(
      fields_split = str_split(fields, "\\|"),
      ALLELE = purrr::map_chr(fields_split, ~ .x[1]),
      SYMBOL = purrr::map_chr(fields_split, ~ .x[2]),
      splice_sites = purrr::map(fields_split, tail, -2)
    ) %>%
    select(-fields, -fields_split) %>%

  # extract and filter individual effects
    tidyr::unnest(splice_sites) %>%
    dplyr::filter(splice_sites != "") %>%
    tidyr::separate(splice_sites, into = c("effect", "score", "pos_rel"), sep = ":") %>%
    mutate(
      score = as.numeric(score),
      pos_rel = as.integer(pos_rel),
      effect = as.factor(effect)
    )

  fix_df %>%
    dplyr::left_join(splice_df, by = "Key")

}
