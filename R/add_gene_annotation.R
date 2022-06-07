#' Annotate splice junctions with gene names.
#'
#' @param df  A data.frame with splice junctions in rows and at least the columns:
#'   -  `junc_id` hash id for a context sequence
#' @param txdb  knownGene Track object such as as `TxDb.Hsapiens.UCSC.hg19.knownGene`
#'
#' @param annotation_db  Genome wide annotation file such as `org.Hs.eg`
#'
#' @return A data.frame as the input with the additional column(s):
#' - `gene_name`: Name of the input sequence.
#'
#' @examples
#'
#'
#'@import dplyr
#'
#'@export
add_gene_name <- function(df, txdb = NULL, annotation_db = NULL){

  stopifnot("junc_id" %in% names(df))

  if(is.null(txdb)){
    message("INFO: Use default txdb known genes from BSgenome.Hsapiens.UCSC.hg19")
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  }

  if(is.null(annotation_db)){
    suppressMessages(require(org.Hs.eg.db))
    message("INFO: Use default annotation from org.Hs.eg")
    annotation_db <- "org.Hs.eg"
  }

  df_gr <- junc_to_gr(unique(df$junc_id))

  # get gene id
  df_id <- suppressMessages(annotate_gene_ids(df_gr, txdb))
  # get gene name
  df_id_name <- df_id %>%
    rowwise()%>%
    mutate(gene_name = ifelse(!purrr::is_empty(gene_id) , toString(annotate::getSYMBOL(na.omit(gene_id), data=annotation_db)), NA)) %>%
    mutate(junc_id = generate_junction_id(seqnames, start, end, strand)) %>%
    dplyr::select(gene_name, junc_id)


  # map gene names on input data
  df_annotated <- df %>%
    left_join(df_id_name, by = "junc_id")


  return(df_annotated)
}


#' Annotate coordinate intervals with gene ids. Modified from Luca Zammator (2019) from http://www.bioinsteps.com/2019/08/annotate-genomic-coordinates-using-r.html
#'
#' @param intervals  A granges object of splice junctions.
#' @param txdb  A txdb objext with known gene annotation such as `TxDb.Hsapiens.UCSC.hg19.knownGene`
#'
#' @return A data.frame with junction range and gene_id
#'
#' @examples
#'
#'
#'@import dplyr
#'
annotate_gene_ids <- function(intervals, txdb){

    stopifnot(is(intervals, "GRanges"), is(txdb, "TxDb"))

    anno <- suppressMessages(GenomicFeatures::genes(txdb))
    olaps <- GenomicRanges::findOverlaps(intervals, anno)
    mcols(olaps)$gene_id <- anno$gene_id[subjectHits(olaps)]
    intervals_factor <- factor(queryHits(olaps), levels = seq_len(queryLength(olaps)))
    intervals$gene_id <-suppressMessages(splitAsList(mcols(olaps)$gene_id, intervals_factor))
    intervals <- as_tibble(intervals)
    return(intervals)

  }
