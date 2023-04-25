
#' Formats spliceAI output and filter for predicted effects
#'
#' Reformat the data for each annotated effect per row, filters effects
#' to have a probability score not NA and score > 0, and removes gene symbol from
#' data to make non-redundant output.
#'
#' @param spliceai_variants [tibble][tibble::tibble-package] with parsed
#' spliceAI mutations from \code{\link{parse_spliceai}}
#' @param gene_table optional [tibble][tibble::tibble-package] with the columns:
#' - `gene_id`: ENSEMBL gene id
#' - `gene_name`: gene symbol
#' If `gene_table` is provided, the formatted data contains a column with the `gene_id`.
#'
#' @return A [tibble][tibble::tibble-package] with splicing effects per row
#'
#' @examples
#' spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")
#' df <- parse_spliceai(spliceai_file)
#' format_spliceai(df)
#'
#' @seealso \code{\link{parse_spliceai}}, \code{\link{annotate_mut_effect}}
#' @export
format_spliceai <- function(spliceai_variants, gene_table = NULL){

  # get all splicing affects for each variant in rows
  spliceai_variants %>%
    pivot_longer(
      cols = DS_AG:DP_DL,
      names_to = c(".value", "effect"),
      names_pattern = "(DS|DP)_(\\w*)",
    ) %>%
    dplyr::rename(score = DS, pos_rel = DP) %>%
    mutate(effect = as.factor(effect)) %>%

    # filter effects without probability score given
    filter(!is.na(score) & score > 0) %>%

    # add unique IDs for mutation
    mutate(
      mut_id = str_c(CHROM, POS, REF, ALT, sep = "_"),
      chr = CHROM,
      pos = as.integer(POS) + pos_rel
    ) %>%

    # if gene_table, annotate gene_id
    # keep only relevant columns
    {
      if (!is.null(gene_table))
        left_join(., gene_table, by = c("SYMBOL" = "gene_name")) %>%
        dplyr::select(mut_id, effect, score, chr, pos_rel, pos, gene_id)
      else
        dplyr::select(. , mut_id, effect, score, chr, pos_rel, pos)
    } %>%

    dplyr::distinct()
}
