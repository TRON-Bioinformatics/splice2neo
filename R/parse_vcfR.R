#' Unified parsing of VCF files using vcfR.
#'
#' @param vcf_file Path to a single VCF.
#' @param info_fields The INFO fields to be extracted into the 'info' df of the
#'  output list.
#' @return A list with two elements: 'info', the info data of the VCF as a
#'  [tibble][tibble::tibble-package] and 'fix', a
#'  [tibble][tibble::tibble-package] of the standard variant fields of the VCF.
#' @keywords internal
parse_vcf <- function(vcf_file, info_fields) {
  vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)

  info_df <- vcf %>%
    vcfR::extract_info_tidy(
      info_fields = info_fields
    )

  fix_df <- vcfR::getFIX(vcf) %>%
    {
      # When only EXACTLY one variant is found `getFIX` will return a
      # character vector.
      if (nrow(vcf) == 1) {
        tibble::as_tibble_row(.)
      } else {
        tibble::as_tibble(.)
      }
    } %>%
    mutate(Key = row_number())

  list(
    "info" = info_df,
    "fix" = fix_df
  )
}
