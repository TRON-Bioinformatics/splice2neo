
#' Annotate if there is an exon within an intron
#'
#' @param junc_id Junction id
#' @param tx_id Transcript id
#' @param patient_id Patient id. This is only relevant if a dataset from multiple patients is annotated.
#' @param transcripts a GRangesList with transcripts defined as GRanges of exons
#'   created by `GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)`.
#'
#' @return Boolean indicating if the given IR junctions is exon-free. Will return NA for non-intron rententions.
#' @examples
#'
#' library(dplyr)
#'
#' toy_tx <- GenomicRanges::GRangesList(list(
#'tx1 = GenomicRanges::GRanges(
#'  c("1", "1", "1"),
#' IRanges::IRanges(
#'    c(2, 15, 27),
#'  c(8, 23, 35),
#'),
#'strand = c("+", "+", "+")
#'),
#'tx1a = GenomicRanges::GRanges(
#'  c("1" ),
#'IRanges::IRanges(
#'  c(3),
#'  c(10),
#'),
#'strand = c("+")
#')))
#'
#'df <- tibble(
#'  junc_id = c("1:8-9:+", "1:14-15:+"),
#'  tx_id = c("tx1", "tx1")
#')
#'df1 <- df %>%
#'  mutate(exon_free = exon_in_intron(junc_id = junc_id, tx_id = tx_id, transcripts = toy_tx))
#'
#' @export
exon_in_intron <- function(junc_id, tx_id, transcripts, patient_id = NULL){

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
    mutate(junc_tx_pat_id = paste0(junc_id, "_", tx_id, "_", patient_id)) %>%
    separate(junc_id,
             into = c("chrom", "start_end" , "strand"),
             sep = ":") %>%
    separate(start_end,
             into = c("junc_start", "junc_end"),
             sep = "-") %>%
    mutate(junc_start = as.integer(junc_start)) %>%
    mutate(junc_end = as.integer(junc_end)) %>%
    mutate(other_position  = as.integer(other_position)) %>%
    mutate(
      new_start = case_when(
        !start_on_exon &
          strand == "+" | start_on_exon & strand == "-" ~ other_position,
        start_on_exon &
          strand == "+"  | !start_on_exon & strand == "-" ~ junc_end
      )
    )

  intron_ranges <- intron_ranges %>%
    mutate(
      new_end = case_when(
        !start_on_exon &
          strand == "+" |  start_on_exon & strand == "-" ~ junc_start,
        start_on_exon &
          strand == "+" | !start_on_exon & strand == "-" ~ other_position,
      )
    ) %>%
    mutate(interval_range = ifelse(new_end != -1, generate_junction_id(chrom, new_start, new_end, strand), NA))

  # get junctions as GRanges object
  jx_df <- intron_ranges %>%
    filter(!is.na(interval_range))


  if(nrow(jx_df) > 0){
    jx <- junc_to_gr(jx_df$interval_range)
  }else{
    jx <- GenomicRanges::GRanges("chr1", IRanges::IRanges(0, 0), "+")
  }


  jx_df <- suppressWarnings(
    jx_df %>%
      mutate(interval_exon_overlap = IRanges::overlapsAny(jx, transcripts , ignore.strand = TRUE)) %>%
      dplyr::select(interval_exon_overlap, junc_tx_pat_id)
  )

  intron_ranges <- intron_ranges %>%
    left_join(jx_df, by = "junc_tx_pat_id")

  # returns TRUE if no exon of other transcript in the given intron region
  # returns FALSE if there is an exon of another transcript
  # returns NA if no intron retention or if the predicted event is at the end or start of a transcript
  exon_free <- case_when(
    abs(intron_ranges$junc_end - intron_ranges$junc_start) != 1 ~ NA,
    intron_ranges$interval_exon_overlap ~ FALSE,
    !intron_ranges$interval_exon_overlap ~TRUE
  )

  return(exon_free)

}
