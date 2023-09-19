# StringTie FUNCTIONS -----------------------------------------------------

#' Imports StringTie GTF file
#'
#' @param file.name Path to StringTie GTF
#'
#' @return A tibble with one splice junction per row and columns
#'  chromosome, junction_start, junction_end, strand, junc_id
#'
#' @import readr stringr
#' @export
import_stringtie_gtf <- function(file.name) {
  if(!file.exists(file.name)){
    stop("StringTie GTF file is missing...")
  }
  message('Importing StringTie gtf...')
  # Import GTF as a GRanges object
  stringtie_gtf <- rtracklayer::import(file.name)
  # Create an auxillary table that holds the exon counts for each
  # transcript/transfrags. This is required to remove single exon
  # fragments that have no strand information and therefore make
  # a generation of a TxDb impossible
  aux_tib <- as_tibble(stringtie_gtf) %>%
    count(gene_id, transcript_id, type, name="typeCount") %>%
    filter(type == "exon") %>%
    filter(typeCount >= 2)
  # Remove single exon transfrags/transcripts
  stringtie_gtf <- subset(stringtie_gtf,
    gene_id %in% aux_tib$gene_id & transcript_id %in% aux_tib$transcript_id)
  # Perform sort of tibble elements. StringTie uses a different way of
  # sorting the GTF file. While GENCODE GTF files sort all the features
  # from upstream to downstream for both strands (ascending on + and descending on -),
  # StringTie sorts features with ascending coordinates regardless of the strand.
  # This sort is required for splice2neo::canonical_junctions
  stringtie_tibble <- as_tibble(stringtie_gtf) %>%
    dplyr::group_by(gene_id, transcript_id) %>%
    dplyr::group_modify(
      ~if(all(.x$strand == '-')){
        dplyr::arrange(.x, desc(start), .by_group=TRUE)
      } else {
        dplyr::arrange(.x, start, .by_group=TRUE)
      }
    ) %>%
    ungroup()
  stringtie_junc <-
    GenomicRanges::makeGRangesListFromDataFrame(
      stringtie_tibble, split.field="transcript_id")
  stringtie_junc <- splice2neo:::canonical_junctions(stringtie_junc)
  fields <- stringr::str_match(stringtie_junc, "(.*):(\\d+)-(\\d+):([*+-])")

  junc_tibble <- tibble(chr = fields[,2],
                        start = fields[,3] %>% as.integer(),
                        end = fields[,4] %>% as.integer(),
                        strand = fields[,5])

  colnames(junc_tibble) <-
    c(
      "chromosome",
      "junction_start",
      "junction_end",
      "strand",
      "junc_id"
    )
  return(junc_tibble)
}

#' Transforms StringTie intermediate table into standardized junction format
#'
#' @param tib StringTie junctions as tibble
#'
#' @return A tibble in standardized junction format
#'
#'
#' @import dplyr
#' @export
stringtie_transform_format <- function(tib) {
  tib <- tib %>%
    mutate(
      class = NA,
      AS_event_ID = NA,
      Gene = NA,
      junc_id = generate_junction_id(chromosome, junction_start, junction_end, strand)
    ) %>%
    sort_columns()
  return(tib)
}

#' Imports StringTie assembled transcripts and transforms the raw output
#' into standardized junction output format
#'
#'  - GitHub: https://github.com/gpertea/stringtie
#'  - Paper: https://doi.org/10.1038/nbt.3122
#'
#' @param gtf.file The path to StringTie GTF file
#'
#' @return A tibble in standardized junction format
#'
#' @export
stringtie_transform <- function(gtf.file) {
  dat <- gtf.file %>%
    import_stringtie_gtf() %>%
    stringtie_transform_format()
  return(dat)
}

