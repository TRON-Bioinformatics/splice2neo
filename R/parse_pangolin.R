
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
  # Format: gene|pos:largest_increase|pos:largest_decrease|,gene|pos:largest_increase|pos:largest_decrease|
  # if pangolin option -s was used, there can be multiple increase/decrease scores!

  col_names <- c("gene_id",
                 "pos_largest_increase","pos_largest_decrease")

  vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)

  # get fixed fields for variants
  fix_df <- vcfR::getFIX(vcf)

  # required for vcf with only one mutation
  if(is.null(nrow(fix_df))){
    fix_df <- fix_df %>%
      t()
  }

  fix_df <- fix_df %>%
    tibble::as_tibble() %>%
    mutate(
      Key = row_number()
    )

  # get annotations
  pangolin_df <- vcf %>%
    vcfR::extract_info_tidy(info_fields = c("Pangolin")) %>%
    mutate(
      annot_gene = Pangolin %>%

        # remove warning massage
        str_replace_all("Warnings:NoAnnotatedSitesToMaskForThisGene", "") %>%
        str_replace_all("Warnings:", "") %>%

        # extract all gene-wise annotations per variant
        str_split(",")
    ) %>%
    select(-Pangolin) %>%
    unnest(annot_gene) %>%
    # pares individual annotation values
    mutate(
      gene_id = gsub("\\|.*", "", annot_gene),
      # extract all annotations per variant - gene
      annot = annot_gene %>%
        str_extract_all("-?\\d+:-?\\d\\.+\\d+")
    ) %>%
    unnest(annot) %>%
    mutate(
      pos_rel = annot %>%
        str_replace("(-?\\d+):-?\\d\\.+\\d+", "\\1") %>%
        as.integer(),
      pangolin_score = annot %>%
        str_replace("-?\\d+:(-?\\d\\.+\\d+)", "\\1") %>%
        as.numeric()
    ) %>%
    select(-annot, -annot_gene)


  fix_df %>%
    left_join(pangolin_df, by = "Key")

}
