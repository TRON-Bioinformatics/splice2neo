
#' Annotate splice variants from spliceAI with resulting junctions
#'
#' The mutation effects donor loss (DL), donor gain (DG), acceptor loss (AL),
#' and acceptor gain (AG) are annotated with respect to the provided transcripts
#' to build resulting splice junctions.
#' Donor loss (DG) and acceptor loss (AD) are only considered when overlapping
#' with an exon-intron junction in the provided transcripts.
#'
#' @param var_df a data.frame wit variants and (at least) the following columns:
#'   - `CHROM`
#'   - `POS`
#'   - `REF`
#'   - `ALT`
#'   - `change` an change effect class from *SpliceAI*. One of `DL`, `DG`, `AL`, `AG`.
#'   - `pos_rel` affected position relative to `POS`.
#'
#' @param transcripts a GRangesList with transcripts defined as GRanges of exons
#'   created by `GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)`.
#' @param transcripts_gr a GRanges object with transcript created by
#'   `GenomicFeatures::transcripts(txdb)`
#'
#' @return A data.frame with with additional rows and columns including the
#' splice junction in the column `junc_id`.
#'
#' @examples
#' spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")
#' df_raw <- parse_spliceai(spliceai_file)
#' df <- format_spliceai(df_raw)
#'
#' annotate_spliceai_junction(df, toy_transcripts, toy_transcripts_gr)
#'
#'
#' @export
annotate_spliceai_junction <- function(var_df, transcripts, transcripts_gr){

  var_df <- var_df %>%
    mutate(
      mut_id = str_c(CHROM, POS, REF, ALT, sep = "_"),
      mut_effect_id = str_c(mut_id, "_", row_number())
    )

  var_gr <- GenomicRanges::GRanges(str_c(
    var_df$CHROM,
    ":",
    as.integer(var_df$POS) + var_df$pos_rel
    ),
    # var_effect_id = var_df$var_effect_id
  )
  names(var_gr) <- var_df$mut_effect_id

  message("INFO: calculate coordinates of upstream and downstream exons...")
  # get all possible junctions of by start and end coordinates of upsteam and downstream exons
  next_junc_df <- next_junctions(var_gr, transcripts, transcripts_gr)

  message("INFO: calculate junction coordinates from predicted effect...")

  junc_df <- var_df %>%

    # add next exon coordinates of next exons
    left_join(next_junc_df, by = "mut_effect_id") %>%

    # Filter out donor loss and acceptor loss which is not on exon-intron boundaries
    filter(
      change != "DL" | at_end,
      change != "AL" | at_start
    ) %>%

    # add rules
    mutate(
      change = as.character(change),
      pos = as.integer(POS) + pos_rel
    ) %>%
    left_join(
      change_to_junction_rules,
      by = c("change")
    ) %>%

    # apply rules
    rowwise() %>%

    # evaluate rule to get coordinates of junctions
    mutate(
      strand_offset = ifelse(tx_strand == "-", -1, 1),
      coord_1 = eval(parse(text = rule_left)),
      coord_2 = eval(parse(text = rule_right)),
      left = ifelse(coord_1 <= coord_2, coord_1, coord_2),
      right = ifelse(coord_1 <= coord_2, coord_2, coord_1),
    ) %>%

    # remove predicted effects with missing values
    filter(!is.na(left) & !is.na(right) & !is.na(tx_strand) & !is.na(CHROM)) %>%
    # remove predicted effects outside of transcript range
    filter((!is.na(downstream_start) & !is.na(downstream_end)) | (!is.na(upstream_start) & !is.na(upstream_end)))%>%

    # add junction IDs
    mutate(
      junc_id = generate_junction_id(CHROM, left, right, tx_strand),
      tx_junc_id = str_c(tx_id, junc_id, sep = "_"),
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
#'  upstream_end, downstream_start, downstream_end, at_start, at_end.
#'
next_junctions <- function(var_gr, transcripts, transcripts_gr){

  # get mut_id
  var_names <- names(var_gr)

  message("INFO: get overlapping transcripts..." )
  # get overlapping transcripts
  hits <- suppressWarnings(GenomicRanges::findOverlaps(var_gr, transcripts_gr))

  message("INFO: Build df ..." )
  var_to_transcript <- hits %>%
    as.data.frame() %>%
    as_tibble() %>%
    dplyr::rename(
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
      mut_effect_id = var_names[var_nr],

      # tx annotation
      tx_chr = as.character(GenomeInfoDb::seqnames(transcripts_gr)[tx_nr]),
      tx_strand = as.character(BiocGenerics::strand(transcripts_gr)[tx_nr]),
      tx_id = if(!is.null(names(transcripts_gr))){names(transcripts_gr)[tx_nr]} else{as.character(transcripts_gr$tx_name)[tx_nr]},

      # extract all exon start and end coordinates of transcripts as GRanges
      starts_gr = map(tx_gr, GenomicRanges::resize, width = 1, fix = "start"),
      ends_gr = map(tx_gr, GenomicRanges::resize, width = 1, fix = "end"),

      # get the closest upstream and downstream start and end positions
      upstream_start_idx = map2_int(var_gr, starts_gr, GenomicRanges::follow),
      downstream_start_idx = map2_int(var_gr, starts_gr, GenomicRanges::precede),
      upstream_end_idx = map2_int(var_gr, ends_gr, GenomicRanges::follow),
      downstream_end_idx = map2_int(var_gr, ends_gr, GenomicRanges::precede),

      upstream_start = map2_int(starts_gr, upstream_start_idx, ~ifelse(!is.na(.y), BiocGenerics::start(.x[.y]), NA)),
      downstream_start = map2_int(starts_gr, downstream_start_idx, ~ifelse(!is.na(.y), BiocGenerics::start(.x[.y]), NA)),
      upstream_end = map2_int(ends_gr, upstream_end_idx, ~ifelse(!is.na(.y), BiocGenerics::start(.x[.y]), NA)),
      downstream_end = map2_int(ends_gr, downstream_end_idx, ~ifelse(!is.na(.y), BiocGenerics::start(.x[.y]), NA)),

      # calculate if effect position overlaps with exon start or exon end
      at_start = suppressWarnings(map2_lgl(var_gr, starts_gr, IRanges::overlapsAny)),
      at_end = suppressWarnings(map2_lgl(var_gr, ends_gr, IRanges::overlapsAny)),

    ) %>%

    dplyr::select(mut_effect_id, var_nr, tx_chr, tx_id, tx_strand, upstream_start,
           upstream_end, downstream_start, downstream_end, at_start, at_end)

}

#' Rules on how a splicing affecting variant creates a junction
change_to_junction_rules <- tribble(
  ~change, ~class,             ~rule_left,        ~rule_right,
  "DL",    "intron retention", "pos",             "pos + strand_offset",
  "DL",    "exon skipping",    "upstream_end",    "downstream_start",

  "DG",    "alternative 5prim", "pos",             "downstream_start",

  "AL",    "intron retention",  "pos - strand_offset",         "pos",
  "AL",    "exon skipping",     "upstream_end",    "downstream_start",

  "AG",    "alternative 3prim",  "upstream_end",     "pos",
)
