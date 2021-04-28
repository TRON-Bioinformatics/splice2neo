
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

}
