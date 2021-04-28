
#' Formats spliceAI output and filter for predicted effects
#'
#' Reformats the data for each annotated effect per row and filters effects
#' to have a probability not NA and score > 0.
#'
#' @param spliceai_variants [tibble][tibble::tibble-package] with parsed
#' spliceAI mutations from \code{\link{parse_spliceai}}
#'
#' @return A [tibble][tibble::tibble-package] with spliceing effects per row
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
    rename(prob = DS, pos_rel = DP) %>%
    mutate(change = as.factor(change)) %>%

    # filter effects without probability given
    filter(!is.na(prob) & prob > 0)
}
