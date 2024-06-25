
#' Formats CI-SpliceAI with -t flag output and filter for predicted effects
#'
#' Reformat the data for each annotated effect per row, filters effects
#' to have a probability score not NA and score > 0, and removes gene symbol from
#' data to make non-redundant output.
#'
#' @param cispliceai_variants [tibble][tibble::tibble-package] with parsed
#' CI-SpliceAI mutations from \code{\link{parse_spliceai}}
#' @param transcripts_gr *Optionally*, A GRanges object with transcript ranges created by
#'   `GenomicFeatures::transcripts(txdb)`can be provided which will allow to annotate full gene ids
#'   to the formatted table and to consider only transcripts related to the annotated gene in `annotate_mut_effect`.
#'   This parameter is optionally.
#'
#' @return A [tibble][tibble::tibble-package] with splicing effects per row
#'
#' @examples
#' cispliceai_file <- system.file("extdata", "cispliceai_thresh_output.vcf", package = "splice2neo")
#' df <- parse_cispliceai_thresh(cispliceai_file)
#' format_cispliceai_thresh(df)
#'
#' @seealso \code{\link{parse_cispliceai_thresh}}, \code{\link{annotate_mut_effect}}
#' @export
format_cispliceai_thresh <- function(cispliceai_variants, transcripts_gr = NULL){

  #format columns
  cispliceai_variants <- cispliceai_variants %>%
    mutate(score = as.numeric(score),
           pos_rel = as.integer(pos_rel),
           effect = as.factor(effect)) %>%

    # filter effects without probability score given
    filter(!is.na(score) & score > 0) %>%

    # add unique IDs for mutation
    mutate(
      mut_id = str_c(CHROM, POS, REF, ALT, sep = "_"),
      chr = CHROM,
      pos = as.integer(POS) + pos_rel
    )

  if(!is.null(transcripts_gr)){

    # CI-SpliceAI does not return the version part of the gene id if gene table from CI-SPliceAi is used as annotation
    # we need to map them based on the transcripts_gr object
    # if gene_id should be kept
    gene_table <- tibble::tibble(gene_id = unique(unlist(transcripts_gr@elementMetadata$gene_id))) %>%
      dplyr::mutate(gene_id_short = gsub("\\..*", "", gene_id))

    formated_variants <- cispliceai_variants %>%
      dplyr::left_join(gene_table, by = c("SYMBOL" = "gene_id_short")) %>%
      dplyr::select(mut_id, effect, score, chr, pos_rel, pos, gene_id)

  } else{

    formated_variants <- cispliceai_variants %>%
      dplyr::select(mut_id, effect, score, chr, pos_rel, pos)

  }

  formated_variants %>%
    dplyr::distinct()
}
