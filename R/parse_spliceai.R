
#' Parse VCF output file from spliceAI as data.frame
#'
#' requires the vcfR package.
#'
#' @param vcf_file path to a single VCF file, the output of spliceAI
#' @return a data.frame like object with one row per variant and colums for
#' annotation
#'
#' @examples
#'
#' spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")
#' parse_spliceai(spliceai_file)
#'
#' @import dplyr stringr tidyr
#'
#' @export
parse_spliceai <- function(vcf_file){

  # Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL

  col_names <- c("ALLELE","SYMBOL","DS_AG",
                 "DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL")

  vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)

  # get fixed fields for variants
  fix_df <-  vcfR::getFIX(vcf) %>%
    tibble::as_tibble() %>%
    mutate(
      Key = row_number()
    )

  # spliceai_df <-
  splice_df <- vcf %>%
    vcfR::extract_info_tidy(info_fields = c("SpliceAI")) %>%
    mutate(
      fields = SpliceAI %>% stringr::str_split(",")
    ) %>%
    select(-SpliceAI) %>%
    unnest(fields) %>%
    separate(fields, into = col_names, sep = "\\|") %>%
    mutate_at(vars(starts_with("DS_")), as.numeric) %>%
    mutate_at(vars(starts_with("DP_")), as.integer)

  fix_df %>%
    left_join(splice_df, by = "Key")

}
