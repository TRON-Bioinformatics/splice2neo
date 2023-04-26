
#' Annotate splice variants effects with resulting junctions
#'
#' The mutation effects donor loss (DL), donor gain (DG), acceptor loss (AL),
#' and acceptor gain (AG) are annotated with respect to the provided transcripts
#' to build resulting splice junctions.
#' Donor loss (DG) and acceptor loss (AD) are only considered when overlapping
#' with an exon-intron junction in the provided transcripts.
#'
#' @param effect_df a data.frame with variant effects on splicing per row. It
#' should have at least the following columns:
#'   - `chr` chromosome
#'   - `pos` absolute position of the effect
#'   - `effect` an splicing effect class, one of `DL`, `DG`, `AL`, `AG`.
#'   optional column:
#'   - `gene_id` ENSEMBL gene id. Is required if `gene_mapping` is set to TRUE
#'
#' @param transcripts a GRangesList with transcripts defined as GRanges of exons
#'   created by `GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)`.
#' @param transcripts_gr a GRanges object with transcript created by
#'   `GenomicFeatures::transcripts(txdb, columns = c("gene_id", "tx_id", "tx_name"))`
#' @param gene_mapping Indicator whether only transcripts related to the gene
#' provided in the `gene_id` column in `effect_df` are considered.
#' If `gene_mapping` is FALSE, all potentially affected  transcripts covering
#' the relevant genomic position are considered.
#' If `gene_mapping` is TRUE, potentially affected transcripts from the gene
#' provided in `effect_df` that cover the relevant genomic positions are considered.
#'
#' @return A data.frame with with additional rows and columns including the
#' splice junction in the column `junc_id`.
#'
#' @examples
#' spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")
#' df_raw <- parse_spliceai(spliceai_file)
#' effect_df <- format_spliceai(df_raw)
#'
#' annotate_mut_effect(effect_df, toy_transcripts, toy_transcripts_gr)
#'
#'
#' @export
annotate_mut_effect <- function(effect_df, transcripts, transcripts_gr, gene_mapping = FALSE){

  effect_df <- effect_df  %>%
    mutate(
      effect_index = str_c("EIDX_", row_number())
    )

  if(nrow(effect_df) != 0){
    var_gr <- GenomicRanges::GRanges(
      str_c(
        effect_df$chr,
        ":",
        as.character(effect_df$pos)
    ))
    names(var_gr) <- effect_df$effect_index

    message("INFO: calculate coordinates of upstream and downstream exons...")
    # get all possible junctions of by start and end coordinates of upsteam and downstream exons
    next_junc_df <- next_junctions(var_gr, transcripts, transcripts_gr)

    message("INFO: calculate junction coordinates from predicted effect...")

    junc_df <- effect_df %>%

      # add next exon coordinates of next exons
      left_join(next_junc_df, by = "effect_index") %>%

      # Filter out donor loss and acceptor loss which is not on exon-intron boundaries
      filter(
        effect != "DL" | at_end,
        effect != "AL" | at_start
      ) %>%

      # add rules
      mutate(
        effect = as.character(effect),
        # pos = as.integer(POS) + pos_rel
      ) %>%
      left_join(
        effect_to_junction_rules,
        by = c("effect")
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
      filter(!is.na(left) & !is.na(right) & !is.na(tx_strand) & !is.na(chr)) %>%
      # remove predicted effects outside of transcript range
      filter((!is.na(downstream_start) & !is.na(downstream_end)) | (!is.na(upstream_start) & !is.na(upstream_end)))%>%

      # add junction IDs
      mutate(
        junc_id = generate_junction_id(chr, left, right, tx_strand),
        tx_junc_id  = str_c(tx_id, junc_id, sep = "_"),
      ) %>%
      ungroup()


    if(gene_mapping){

      if(!"gene_id" %in% names(junc_df)) {
        stop("The column gene_id is required in effect_df if gene mapping should be applied")
      }

      # gene and transcripts
      gene_transcript_mapping <-
        tibble::tibble(
          gene_id = unlist(transcripts_gr@elementMetadata$gene_id),
          tx_name = transcripts_gr@elementMetadata$tx_name
        )

      # consider only transcripts that relate to gene that is provided in the
      # output of mutation effect prediction tools
      junc_df <- junc_df %>%
        dplyr::left_join(
          gene_transcript_mapping,
          by = c("tx_id" = "tx_name"),
          suffix = c("", "_mapped")
        ) %>%
        group_by(mut_id, junc_id, effect) %>%
        dplyr::filter(gene_id == gene_id_mapped) %>%
        dplyr::select(-gene_id_mapped) %>%
        ungroup()

    }

    message("INFO: Evaluation of rules done.")
  } else{
    message("WARNING: There are no mutations with predicted splice effect by SpliceAI")
    junc_df <- effect_df %>%
      tibble::add_column(
        "var_nr" = NA ,
        "tx_chr"= NA,
        "tx_id" = NA,
        "tx_strand" =NA,
        "upstream_start"= NA,
        "upstream_end" = NA,
        "downstream_start"= NA,
        "downstream_end"= NA,
        "at_start"= NA,
        "at_end"= NA,
        "event_type"= NA,
        "rule_left"= NA,
        "rule_right"= NA,
        "strand_offset"= NA,
        "coord_1"= NA ,
        "coord_2"= NA,
        "left"= NA,
        "right" = NA,
        "junc_id"= NA,
        "tx_junc_id"= NA
      )
  }

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
      var_gr = purrr::map(var_nr, ~var_gr[.x]),
      tx_gr = transcripts[tx_nr] %>% as.list()
    )

  message("INFO: find overlapping exons..." )

  var_to_transcript <- var_to_transcript %>%
    mutate(

      # var annot
      effect_index = var_names[var_nr],

      # tx annotation
      tx_chr = as.character(GenomeInfoDb::seqnames(transcripts_gr)[tx_nr]),
      tx_strand = as.character(BiocGenerics::strand(transcripts_gr)[tx_nr]),
      tx_id = if(!is.null(names(transcripts_gr))){names(transcripts_gr)[tx_nr]} else{as.character(transcripts_gr$tx_name)[tx_nr]},

      # extract all exon start and end coordinates of transcripts as GRanges
      starts_gr = purrr::map(tx_gr, GenomicRanges::resize, width = 1, fix = "start"),
      ends_gr = purrr::map(tx_gr, GenomicRanges::resize, width = 1, fix = "end"),

      # get the closest upstream and downstream start and end positions
      upstream_start_idx =purrr:: map2_int(var_gr, starts_gr, GenomicRanges::follow),
      downstream_start_idx = purrr::map2_int(var_gr, starts_gr, GenomicRanges::precede),
      upstream_end_idx = purrr::map2_int(var_gr, ends_gr, GenomicRanges::follow),
      downstream_end_idx = purrr::map2_int(var_gr, ends_gr, GenomicRanges::precede),

      upstream_start = purrr::map2_int(starts_gr, upstream_start_idx, ~ifelse(!is.na(.y), BiocGenerics::start(.x[.y]), NA)),
      downstream_start = purrr::map2_int(starts_gr, downstream_start_idx, ~ifelse(!is.na(.y), BiocGenerics::start(.x[.y]), NA)),
      upstream_end = purrr::map2_int(ends_gr, upstream_end_idx, ~ifelse(!is.na(.y), BiocGenerics::start(.x[.y]), NA)),
      downstream_end = purrr::map2_int(ends_gr, downstream_end_idx, ~ifelse(!is.na(.y), BiocGenerics::start(.x[.y]), NA)),

      # calculate if effect position overlaps with exon start or exon end
      at_start = suppressWarnings(purrr::map2_lgl(var_gr, starts_gr, IRanges::overlapsAny)),
      at_end = suppressWarnings(purrr::map2_lgl(var_gr, ends_gr, IRanges::overlapsAny)),

    ) %>%

    dplyr::select(effect_index, var_nr, tx_chr, tx_id, tx_strand, upstream_start,
           upstream_end, downstream_start, downstream_end, at_start, at_end)

}

#' Rules on how a splicing affecting variant creates a junction
effect_to_junction_rules <- tribble(
  ~effect, ~event_type,             ~rule_left,        ~rule_right,
  "DL",    "intron retention", "pos",             "pos + strand_offset",
  "DL",    "exon skipping",    "upstream_end",    "downstream_start",

  "DG",    "alternative 5prim", "pos",             "downstream_start",

  "AL",    "intron retention",  "pos - strand_offset",         "pos",
  "AL",    "exon skipping",     "upstream_end",    "downstream_start",

  "AG",    "alternative 3prim",  "upstream_end",     "pos",
)
