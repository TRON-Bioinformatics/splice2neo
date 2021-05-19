
#' Annotate splice variants with resulting junctions
#'
#'
#' @param var_df a data.frame wit variants and (at least) the following columns:
#'   - `CHROM`
#'   - `POS`
#'   - `ALT`
#'   - `ALT`
#'   - `change` an change effect class from *SpliceAI*. One of `DL`, `DG`, `AL`, `AG`.
#'   - `pos_rel` affected position relative to `POS`.
#'
#' @param transcripts a GRangesList with transcripts defined as GRanges of exons
#'   created by `GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)`.
#' @param transcripts_gr a GRanges object with transcript created by
#'   `GenomicFeatures::transcripts(txdb)`
#'
annotate_spliceai_junction <- function(var_df, transcripts, transcripts_gr){

  var_df <- var_df %>%
    mutate(
      mut_id = str_c(CHROM, POS, REF, ALT, sep = "_")
    )

  var_gr <- GenomicRanges::GRanges(str_c(var_df$CHROM, ":",
                          as.integer(var_df$POS) + var_df$pos_rel))
  names(var_gr) <- var_df$mut_id

  message("INFO: calculate coordinates of upstream and downstream exons...")
  # get all possible junctions of by start and end coordinates of upsteam and downstream exons
  next_junc_df <- next_junctions(var_gr, transcripts, transcripts_gr)

  message("INFO: calculate junction coordinates from predicted effect...")

  junc_df <- var_df %>%

    # add next exon coordinates of next exons
    left_join(next_junc_df, by = "mut_id") %>%

    # add rules
    mutate(
      change = as.character(change),
      pos = as.integer(POS) + pos_rel
    ) %>%
    left_join(
      change_to_junction_rules,
      by = c("change", "tx_strand" = "strand")
    ) %>%

    # apply rules
    rowwise() %>%

    # evaluate rule to get coordinates of junctions
    mutate(
      left = eval(parse(text = rule_left)),
      right = eval(parse(text = rule_right)),
    ) %>%

    # add junction IDs
    mutate(
      junc_id = str_c(CHROM, left, right, tx_strand, sep = "_"),
      tx_junc_id = str_c(tx_id, CHROM, left, right, tx_strand, sep = "_"),
    ) %>%
    ungroup()

  message("INFO: Evaluation of rules done.")
  return(junc_df)
}

#' Get genomic coordinates of possible next donor and acceptor sides.
#'
#' @param var_gr a GenomicRanges object with variants (of length one).
#' It is assumed that it is named with a unique ID.
#' Usually this is `<chr>_<pos>_<alt>`.
#'
#' @param transcripts
#'
#' @param transcripts a GRangesList with transcripts defined as GRanges of exons
#'   created by `GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)`.
#' @param transcripts_gr a GRanges object with transcript ranges created by
#'   `GenomicFeatures::transcripts(txdb)`
#'
#' @return A data.frame with possible upstream and downstream exon coordinates for
#' all overlapping transcripts. The data frame contains the following columns:
#'  mut_id, var_nr, tx_chr, tx_id, exon_idx, tx_strand, upstream_start,
#'  upstream_end, downstream_start, downstream_end.
#'
next_junctions <- function(var_gr, transcripts, transcripts_gr){

  # get mut_id
  var_names <- names(var_gr)

  message("INFO: get overlapping transcripts..." )
  # get overlapping transcripts
  hits <- suppressWarnings(GenomicRanges::findOverlaps(var_gr, transcripts))

  message("INFO: Build df ..." )
  var_to_transcript <- hits %>%
    as.data.frame() %>%
    as_tibble() %>%
    rename(
      var_nr = queryHits,
      tx_nr = subjectHits
    )

  message("INFO: build pairs..." )

  var_to_transcript <- var_to_transcript %>%
    mutate(
      #var_gr = map(var_nr, ~var_gr[.x]),
      var_gr = furrr::future_map(var_nr, ~var_gr[.x]),
      # tx_gr = map(tx_nr, ~transcripts[[.x]]),
      tx_gr = transcripts[tx_nr] %>% as.list()
    )

  message("INFO: find overlapping exons..." )

  var_to_transcript <- var_to_transcript %>%
    mutate(

      # var annot
      mut_id = var_names[var_nr],

      # tx annotation
      tx_chr = as.character(GenomeInfoDb::seqnames(transcripts_gr)[tx_nr]),
      tx_strand = as.character(BiocGenerics::strand(transcripts_gr)[tx_nr]),
      tx_id = as.character(transcripts_gr$tx_name)[tx_nr],

      # get upstream donor and downstream acceptor
      exon_idx = furrr::future_map2_int(tx_gr, var_gr,
                                        ~which(suppressWarnings(
                                          IRanges::overlapsAny(.x, .y)
                                          ))),

      # get the next upstream exon (if exists) and take the end coordinate
      upstream_start = furrr::future_map2_int(tx_gr, exon_idx, ~ifelse(.y > 1, IRanges::start(.x)[.y - 1], NA)),
      upstream_end = furrr::future_map2_int(tx_gr, exon_idx, ~ifelse(.y > 1, IRanges::end(.x)[.y - 1], NA)),
      downstream_start = furrr::future_map2_int(tx_gr, exon_idx, ~ifelse(.y < length(.x), IRanges::start(.x)[.y + 1], NA)),
      downstream_end = furrr::future_map2_int(tx_gr, exon_idx, ~ifelse(.y < length(.x), IRanges::end(.x)[.y + 1], NA))
    )  %>%

    select(mut_id, var_nr, tx_chr, tx_id, exon_idx, tx_strand, upstream_start,
           upstream_end, downstream_start, downstream_end)

}

#' Rules on how a splicing affecting variant creates a junction
change_to_junction_rules <- tribble(
  ~change, ~class,             ~strand, ~rule_left,        ~rule_right,
  "DL",    "intron retention", "+",     "pos",             "pos + 1",
  "DL",    "intron retention", "-",     "pos - 1",         "pos",
  "DL",    "exon skipping",    "+",     "upstream_end",    "downstream_start",
  "DL",    "exon skipping",    "-",     "downstream_end",  "upstream_start",

  "DG",    "alternative 5prim", "+",    "pos",             "downstream_start",
  "DG",    "alternative 5prim", "-",    "downstream_end",  "pos",

  "AL",    "intron retention", "+",     "pos - 1",         "pos",
  "AL",    "intron retention", "-",     "pos",             "pos + 1",
  "AL",    "exon skipping",    "+",     "upstream_end",    "downstream_start",
  "AL",    "exon skipping",    "-",     "downstream_end",  "upstream_start",

  "AG",    "alternative 3prim", "+",    "upstream_end",     "pos",
  "AG",    "alternative 3prim", "-",    "pos",              "upstream_start",
)

