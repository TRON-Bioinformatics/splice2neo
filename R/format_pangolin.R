
#' Formats pangolin output and filter for predicted effects
#'
#' Reformat the data for each annotated effect per row. As pangolin does not
#' provide information on donor or acceptor, use both for gain and loss annotations.
#' Effects with score = 0 are filtered out.
#'
#' @param variants [tibble][tibble::tibble-package] with parsed
#' pangolin mutations from \code{\link{parse_pangolin}}
#' @param keep_gene_id Indicator whether the gene_id should be kept in the formatted data.
#'
#' @return A [tibble][tibble::tibble-package] with splicing effects per row
#'
#' pangolin_file <- system.file("extdata", "spliceai_output.pangolin.vcf", package = "splice2neo")
#' df <- parse_pangolin(pangolin_file)
#' format_pangolin(df)
#'
#' @seealso \code{\link{parse_spliceai}}, \code{\link{annotate_mut_effect}}, \code{\link{format_spliceai}}
#' @export
format_pangolin <- function(variants, keep_gene_id = TRUE){

  # get all splicing affects for each variant in rows
  variants %>%
    pivot_longer(
      cols = c("increase_pos", "increase_score", "decrease_pos", "decrease_score"),
      names_sep = "_",
      # names_to = c("effect", "score_pos"),
      names_to = c("effect_direction", ".value"),
      # names_pattern = "(DS|DP)_(\\w*)",
    ) %>%
    dplyr::rename(pos_rel = pos) %>%

    # filter out effects without score or scores <= 0
    filter(!is.na(score) & score > 0) %>%

    # as no annotation on acceptor or donor existst in pangolin, both are considered
    left_join(
      pangolin_effect_translation,
      by = "effect_direction"
    ) %>%

    # add unique IDs for mutation
    mutate(
      mut_id = str_c(CHROM, POS, REF, ALT, sep = "_"),
      chr = CHROM,
      pos = as.integer(POS) + pos_rel
    ) %>%

    # keep only relevant columns
    {
      if (keep_gene_id)
        dplyr::select(. ,mut_id, effect, score, chr, pos_rel, pos, gene_id)
      else
        dplyr::select(. , mut_id, effect, score, chr, pos_rel, pos)
    } %>%

    dplyr::distinct()
}

#' This dataset translates the increase/decrease splicing score from pangolin
#' into a donor gain/loss and acceptor gain/loss effect annotations.
pangolin_effect_translation <- dplyr::tibble(
  effect_direction = c("increase", "increase", "decrease", "decrease"),
  effect = c("DG", "AG", "DL", "AL")
)
