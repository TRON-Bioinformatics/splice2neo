
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
      annot = Pangolin %>%

        # remove warning massage
        str_replace_all("Warnings:NoAnnotatedSitesToMaskForThisGene", "") %>%

        # extract each annotation per variant
        str_extract_all("[^|]*\\|-?\\d+:\\d\\.+\\d+\\|-?\\d+:\\d\\.+\\d+\\|")
    ) %>%
    select(-Pangolin) %>%
    unnest(annot) %>%
    # pares individual annotation values
    mutate(
      gene_id = annot %>%
        str_replace("([^|]*)\\|-?\\d+:\\d\\.+\\d+\\|-?\\d+:\\d\\.+\\d+\\|", "\\1"),
      increase_pos = annot %>%
        str_replace("[^|]*\\|(-?\\d+):\\d\\.+\\d+\\|-?\\d+:\\d\\.+\\d+\\|", "\\1") %>%
        as.integer(),
      increase_score = annot %>%
        str_replace("[^|]*\\|-?\\d+:(\\d\\.+\\d+)\\|-?\\d+:\\d\\.+\\d+\\|", "\\1") %>%
        as.numeric(),
      decrease_pos = annot %>%
        str_replace("[^|]*\\|-?\\d+:\\d\\.+\\d+\\|(-?\\d+):\\d\\.+\\d+\\|", "\\1") %>%
        as.integer(),
      decrease_score = annot %>%
        str_replace("[^|]*\\|-?\\d+:\\d\\.+\\d+\\|-?\\d+:(\\d\\.+\\d+)\\|", "\\1") %>%
        as.numeric(),
    )


  fix_df %>%
    left_join(pangolin_df, by = "Key")

}
