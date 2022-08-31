
#' Annotate if there is an exon within an intron
#'
#' @param df a data.frame wit variants and (at least) the following columns:
#'   - `junc_id`
#'   - `tx_id`
#'
#' @param transcripts a GRangesList with transcripts defined as GRanges of exons
#'   created by `GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)`.
#'
#' @return A data.frame with with an additional row `exon_free`
#'
#'
#'
#' @export
exon_in_intron <- function(junc_id, tx_id, transcripts){



  # df <- tibble(
  #   junc_id = c("chr11:9450577-9450578:+", "chr11:9450176-9450177:+",  "chr6:32150617-32150618:-", "chr6:32150978-32150979:-"),
  #   tx_id = c("ENST00000379719.8_2", "ENST00000379719.8_2", "ENST00000375067.7_2", "ENST00000375067.7_2")
  # )

  # get junctions as GRanges object
  jx <- junc_to_gr(junc_id)
  junc_strand <- str_sub(junc_id,-1,-1)

  # test if junction is an intron retention event
  # TODO: are non IR junctions possible that follow the rule chr:pos-pos+1:strand?
  intron_retention <- ifelse(jx@ranges@width == 2, TRUE, FALSE)

  # get GRanges as subset of transcripts
  tx_lst <- transcripts[tx_id]

  # modify transcripts by appliing the splice junctions
  tx_mod <- modify_tx(tx_lst, jx)

  # only relevant for intron retention events: get the other position of retained
  # interval
  other_genomic_position <- get_intronretention_genomic_alt_pos(tx_lst, tx_mod, jx, intron_retention, junc_strand)
  other_genomic_position$junc_id = junc_id
  other_genomic_position$tx_id = tx_id

  # get interval as id
  intron_ranges <- other_genomic_position %>%
    mutate(junc_tx_id = paste0(junc_id, "_", tx_id)) %>%
    separate(junc_id, into = c("chrom", "start_end" ,"strand"), sep = ":") %>%
    separate(start_end, into = c("junc_start", "junc_end"), sep = "-") %>%
    mutate(junc_start = as.integer(junc_start)) %>%
    mutate(junc_end = as.integer(junc_end)) %>%
    mutate(other_position  = as.integer(other_position )) %>%
    mutate(new_start = case_when(
      !start_on_exon & strand == "+" | start_on_exon & strand == "-" ~ other_position,
      start_on_exon & strand == "+"  | !start_on_exon & strand == "-" ~ junc_end
    )) %>%
    mutate(new_end = case_when(
      !start_on_exon & strand == "+" |  start_on_exon & strand == "-" ~ junc_start,
      start_on_exon & strand == "+" | !start_on_exon & strand == "-" ~ other_position,
    )) %>%
    mutate(interval_range = generate_junction_id(chrom, new_start, new_end, strand))

  # get junctions as GRanges object
  jx_df <- intron_ranges %>%
    filter(!is.na(interval_range))

  jx <- junc_to_gr(jx_df$interval_range)

  jx_df <- jx_df %>%
    mutate(interval_exon_overlap = overlapsAny(jx, transcripts )) %>%
    dplyr::select(interval_exon_overlap, junc_tx_id)

  intron_ranges <- intron_ranges %>%
    left_join(jx_df, by = "junc_tx_id")

  exon_free <- case_when(
    abs(intron_ranges$junc_end - intron_ranges$junc_start) != 1 ~ NA,
    intron_ranges$interval_exon_overlap ~ FALSE,
    !intron_ranges$interval_exon_overlap ~TRUE
  )

  return(exon_free)


}
