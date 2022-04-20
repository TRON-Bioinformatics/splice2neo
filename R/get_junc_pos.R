
#' Get the position of input junction in the transcript sequences
#'
#' @param tx transcripts a \code{\link[GenomicRanges]{GRangesList}} with
#' transcripts defined as GRanges of exons.
#' @param jx splice junctions as GRanges objects
#'
#' @return an integer vector with junction positions
#'
#' @export
get_junc_pos <- function(tx, jx){

  # assume same length of tx and jx
  stopifnot(length(tx) == length(jx))

  # convert to GRangesList
  if(! class(tx) %in% c("GRangesList", "CompressedGRangesList")){
    tx <- GenomicRanges::GRangesList(tx)
  }

  # get individual GRanges objects for start and end position of junction
  j_start <- GenomicRanges::resize(jx, width = 1, fix="start")

  # build common seqlevels
  common_seqlevels <- unique(c(
    GenomeInfoDb::seqlevels(j_start),
    GenomeInfoDb::seqlevels(tx),
    "Z" # add custom seqlevel for zero-size GRangs later
  ))

  GenomeInfoDb::seqlevels(j_start) <- common_seqlevels
  GenomeInfoDb::seqlevels(tx) <- common_seqlevels

  # check if junction overlap the whole range of the transcript ---------------

  # buid full range of transcript
  range_tx <- base::range(tx)

  # check which transcripts have zero length, i.e. consists only of empty ranges
  # these empty ranges get lost when converting to GRanges obejects and therefore
  # will be replaced by an zeor-size toy range object
  empty_ranges <- range_tx %>% GenomicRanges::width() %>% as.integer() %>% is.na()
  if (any(empty_ranges)){
    range_tx[which(empty_ranges)] <-  GenomicRanges::GRangesList(
      GenomicRanges::GRanges("Z", IRanges::IRanges(0, 0), "+")
    )
  }

  range_tx <-  unlist(range_tx)

  on_tx <- IRanges::poverlaps(j_start, range_tx) %>%
    as.logical()

  # in case of empty ranges, on_tx is given as TRUE but should be FALSE
  empty_range <- GenomicRanges::GRanges("Z", IRanges::IRanges(0, 0), "+")
  on_tx <- ifelse(j_start == empty_range, FALSE, on_tx)

  # initialize with NA
  pos_tx <- rep(NA, length(tx))

  # map junction positions on transcript sequence position
  suppressWarnings(
    pos_tx[on_tx] <- GenomicFeatures::pmapToTranscripts(j_start[on_tx], tx[on_tx]) %>%
      BiocGenerics::start()
  )

  return(pos_tx)
}


#' Get the alternative position of input junction from an intron retention event
#'  in the transcript sequences
#'
#' @param tx_lst transcripts a \code{\link[GenomicRanges]{GRangesList}} with
#' transcripts defined as GRanges of exons.
#' @param tx_mod modified transcripts a \code{\link[GenomicRanges]{GRangesList}} with
#' transcripts defined as GRanges of exons.
#' @param jx splice junctions as GRanges objects
#' @param intron_retention logical vector indicating if junctions belong to an intron retention event
#'
#' @return an integer vector with alternative junction positions of intron
#' retention event
#'
#' @export
get_intronretention_alt_pos <- function(tx_lst, tx_mod, jx, intron_retention){

  # assume same length of tx_lst, tx_mod and jx
  stopifnot(length(tx_lst) == length(jx))
  stopifnot(length(tx_mod) == length(jx))
  stopifnot(length(tx_mod) == length(tx_lst))

  # get position of the next exon for intron retentions
  # this is only relevant for intron rententions
  # we need to identify the other boundary of the intron
  jx_lst <- split(jx, 1:length(jx))
  # test if start or and end of junc is on exon
  # get left and right exon indices
  jx_end <- suppressWarnings(mapply(function(j, l) {
    GenomicRanges::findOverlaps(GenomicRanges::resize(j, width = 1, fix="end"), l, select="first")
  } , j = jx_lst, l = as.list(tx_lst), SIMPLIFY = TRUE))
  jx_start <- suppressWarnings(mapply(function(j, l) {
    GenomicRanges::findOverlaps(GenomicRanges::resize(j, width = 1, fix="start"), l, select="first")
  } , j = jx_lst, l = as.list(tx_lst), SIMPLIFY = TRUE))

  start_on_exon <- case_when(is.na(jx_end) & !is.na(jx_start) & intron_retention ~ TRUE,
                             !is.na(jx_end) & is.na(jx_start) & intron_retention ~FALSE
                             )
  empty_range <- GenomicRanges::GRanges("Z", IRanges::IRanges(0, 0), "+")
  # get unknown exon
  ex_range <- mapply(function(j, l, s) {
    if(is.na(s)){ empty_range } else{
      if(s) {ind <- GenomicRanges::precede(j, l)} else if(!s){ind <- GenomicRanges::follow(j, l)}
      if(!is.na(ind)){ l[ind] }else{ empty_range }
    }
  } , j = jx_lst, l = as.list(tx_lst), s = start_on_exon, SIMPLIFY = TRUE)
  ex_range <- suppressWarnings(unlist(as(ex_range, "GRangesList")))

  ex_start <- GenomicRanges::resize(ex_range, width = 1, fix="start")
  ex_end <- GenomicRanges::resize(ex_range, width = 1, fix="end")
  # get position of other boundary of intron
  junc_end_tx = case_when(
    !start_on_exon ~ get_junc_pos(tx_mod, ex_end) ,
    start_on_exon ~ get_junc_pos(tx_mod, ex_start))

  return(junc_end_tx)
}
