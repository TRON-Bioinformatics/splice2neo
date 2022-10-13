
#' Parse VCF output file from pangolin as table
#'
#' requires the vcfR package.
#'
#' @param vcf_file path to a single VCF file, the output of pangolin
#' @return a [tibble][tibble::tibble-package] with one row per variant
#' annotation. Each input variant can have multiple annotations.
#'
#' @examples
#'
#' pangolin_file <- system.file("extdata", "spliceai_output.pangolin.vcf", package = "splice2neo")
#' parse_pangolin(pangolin_file)
#'
#' @import dplyr stringr tidyr
#' @export
parse_pangolin <- function(vcf_file){

  # Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL
  # Format: Pangolin=ENSG00000237298.10_8|33:0.0|-50:0.0|Wa
  # Format: gene|pos:largest_increase|pos:largest_decrease|

  col_names <- c("gene_id",
                 "pos_largest_increase","pos_largest_decrease")

  vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)

  # get fixed fields for variants
  fix_df <-  vcfR::getFIX(vcf) %>%
    tibble::as_tibble() %>%
    mutate(
      Key = row_number()
    )

  # get annotations
  pangolin_df <- vcf %>%
    vcfR::extract_info_tidy(info_fields = c("Pangolin")) %>%
    mutate(
      fields = Pangolin %>%
        str_replace_all("Warnings:NoAnnotatedSitesToMaskForThisGene", "") %>%
        stringr::str_split("\\|")
    ) %>%
    select(-Pangolin) %>%
    unnest(fields) %>%
    filter(!is.na(fields) & fields != "") %>%
    mutate(annot = rep(col_names, nrow(.)/length(col_names))) %>%
    pivot_wider(names_from = annot, values_from = fields)
    # TODO

    separate(fields, into = col_names, sep = "\\|") %>%
    mutate(

      # replace missing data with "." with NA
      across(starts_with("DS_"), na_if, "."),
      across(starts_with("DP_"), na_if, "."),

      # convert scores as numeric and positions as integers
      across(starts_with("DS_"), as.numeric),
      across(starts_with("DP_"), as.integer)
    )

  fix_df %>%
    left_join(splice_df, by = "Key")

}
