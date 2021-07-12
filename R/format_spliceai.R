
#' Formats spliceAI output and filter for predicted effects
#'
#' Reformat the data for each annotated effect per row, filters effects
#' to have a probability not NA and score > 0, and removes gene symbol from
#' data to make non-redundant output.
#'
#' @param spliceai_variants [tibble][tibble::tibble-package] with parsed
#' spliceAI mutations from \code{\link{parse_spliceai}}
#'
#' @return A [tibble][tibble::tibble-package] with splicing effects per row
#'
#' @examples
#' spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")
#' df <- parse_spliceai(spliceai_file)
#' format_spliceai(df)
#'
#' @seealso \code{\link{parse_spliceai}}, \code{\link{annotate_spliceai_junction}}
#' @export
format_spliceai <- function(spliceai_variants){

  # get all splicing affects for each variant in rows
  spliceai_variants %>%
    pivot_longer(
      cols = DS_AG:DP_DL,
      names_to = c(".value", "change"),
      names_pattern = "(DS|DP)_(\\w*)",
    ) %>%
    dplyr::rename(prob = DS, pos_rel = DP) %>%
    mutate(change = as.factor(change)) %>%

    # filter effects without probability given
    filter(!is.na(prob) & prob > 0) %>%

    # remove gene symbol and make non-redundant output
    dplyr::select(-c("SYMBOL")) %>%
    dplyr::distinct()
}
