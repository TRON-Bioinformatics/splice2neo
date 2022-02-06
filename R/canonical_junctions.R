
#' Build canonical junctions from transcripts
#'
#' @param tx a \code{\link[GenomicRanges]{GRangesList}} of reference transcripts
#'
#' @return a character vector of canonical splice junction ids
#'
#'
#' We build all canonical splice junctions that are in the annotated input transcripts.
#' The following lists implementation rules for adjacent canonical exon-exon junction:
#'
#' - strand = +: \eqn{e_i}, \eqn{s_{i+1}}
#' - strand = -: \eqn{e_{i+1}}, \eqn{s_{i}}
#'
#' We also include canonical intron-retention junctions. These are 5' donor or
#' 3' acceptor sites of canonical exon-exon junctions that are not used in all
#' isoforms of the gene. They are located within an exon of other transcripts.
#' Canonical intron-retention junctions are defined by the coordinate of the
#' last exon base and the next base.
#' Therefore, we just need to check whether both bases are included in a single exon.
#'
#'
canonical_junctions <- function(tx){

  # full ranges of all transcript as a single GRanges object
  tx_gr <- base::range(tx)

  exons_gr <- unlist(tx)

  # get data.frame with all transcripts
  tx_df <- dplyr::tibble(
    chr = GenomeInfoDb::seqnames(tx_gr) %>% as.character(),
    strand = BiocGenerics::strand(tx_gr) %>% as.character(),
    exon_start = BiocGenerics::start(tx) %>% as.list(),
    exon_end = BiocGenerics::end(tx) %>% as.list(),
    tx_name = names(tx)
  ) %>%
    unnest(c(exon_start, exon_end))

  # build exon-exon junctions ==================================================
  junc_ee <- tx_df %>%
    dplyr::group_by(tx_name) %>%

    # get left and right junction coordinates as exon boundaries (depending on strand)
    dplyr::mutate(
      left = ifelse(strand == "+", exon_end, dplyr::lead(exon_end)),
      right = ifelse(strand == "+", dplyr::lead(exon_start), exon_start)
    ) %>%
    # filter out entries for first and last exon (without junctions)
    dplyr::filter(!is.na(left), !is.na(right)) %>%
    dplyr::ungroup() %>%
    # - junc_id: `<chr>_<pos1>_<pos2>_<strand>`
    dplyr::mutate(
      junc_id = stringr::str_c(chr, as.integer(left), as.integer(right), strand, sep = "_")
    )

  # build intron-retention junctions ===========================================
  junc_ir <- tx_df %>%

    # get left and right junction coordinates as exon boundaries
    dplyr::mutate(

      start_junc_id = stringr::str_c(chr, exon_start - 1, exon_start, strand, sep = "_"),
      end_junc_id = stringr::str_c(chr, exon_end, exon_end + 1, strand, sep = "_"),

    ) %>%
    # combine into one column
    tidyr::pivot_longer(cols = c("start_junc_id", "end_junc_id"),
                 names_to = "left_right", values_to = "junc_id") %>%
    # distribute junc_id values in to separate columns
    tidyr::separate(junc_id, into = c("chr", "start", "end", "strand"), sep = "_", remove = FALSE)

  # build GRanges object
  junc_ir_gr <- GenomicRanges::GRanges(
    junc_ir$chr,
    IRanges::IRanges(
      as.integer(junc_ir$start),
      as.integer(junc_ir$end),
    ),
    strand = junc_ir$strand,
    junc_id = junc_ir$junc_id
  )

  # get subset of exon start and end junctions that are within an exon
  covered_junc_gr <- junc_ir_gr %>%
    IRanges::subsetByOverlaps(exons_gr, type = "within")

  covered_junc_id <- unique(covered_junc_gr$junc_id)

  # filter to contain only
  junc_ir_covered <- junc_ir %>%
    dplyr::filter(junc_id %in% covered_junc_id) %>%
    dplyr::mutate(
      left = as.integer(start),
      right = as.integer(end)
    ) %>%
    dplyr::select(chr, strand, tx_name, left, right, junc_id)

  # merge with other exon-to-exon junctions ====================================
  junc_df <- dplyr::bind_rows(
    "exon_to_exon" = junc_ee,
    "within_exon" = junc_ir_covered,
    .id = "junction_type"
  )

  # use only unique junctions
  junc_unique <- junc_df %>%
    dplyr::distinct(junc_id, chr, strand, left, right, junction_type) %>%
    dplyr::mutate(in_gencode = TRUE)

  # return only a chracter vector with junctions
  return(unique(junc_df$junc_id))

}
