#' test if junction is canonical junction
#'
#' @param junc_id vector of junction id to test
#' @param ref_junc vector of canonical reference junctions
#' @param exons_gr GRanges of canonical exons
#'
#' @return logical vector of same length as `junc_id` indicating
#' if the input `junc_id` is in the set of reference junctions `ref_junc` or
#' consists of two directly adjacent positions (intron retention) and
#' overlaps (completely within) a canonical exon from `exons_gr`.
#'
is_canonical <- function(junc_id, ref_junc, exons_gr){

  in_ref_junc <- junc_id %in% ref_junc

  # get junction range
  junc_gr <- junc_id_to_gr(junc_id)

  in_exon <- IRanges::overlapsAny(junc_gr, exons_gr, type = "within")
  is_adjacent <- BiocGenerics::width(junc_gr) == 2

  return(in_ref_junc | (in_exon & is_adjacent))

}
