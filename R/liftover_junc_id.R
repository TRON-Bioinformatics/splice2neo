
# TODO:
# - min start, max end for multiple junction liftovers
# - check if both (start and end) coordinates are on same base (length = 1)


# see: https://gitlab.rlp.net/tron/bnt_neoants/neoants/-/blob/dev/workflow/scripts/liftover_junc_id.R?ref_type=heads
# See:  https://gitlab.rlp.net/tron/bnt_neoants/neoants/-/blob/dev/workflow/scripts/check_liftover_sequence.R?ref_type=heads


#' LiftOver junction IDs using the liftOver tool
#'
#' This can be used to liftOver junctions to personalized genome coordinates.
#'
#' @param junc_df a data.frame with at least one column `junc_id` containing junction IDs
#' @param chain_file path to a chain file for the UCSC liftOver tool. See also \code{\link[rtracklayer]{liftOver}}
#' @return a data.frame like the input `junc_df` with the following additional columns:
#'   - `junc_id_lifted_unique` a logical vector indicating if the liftOver was successful and unique (1-to-1 correspondance).
#'   - `junc_id_lifted` a character vector with the lifted junction IDs.
#'   Multiple IDs are separated by `|`.
#'   NA represnt junc_ids that could not be lifted.
#'
#' @examples
#'
#' chain_file = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
#' junc_df <- toy_junc_df
#' liftover_junc_id(junc_df, chain_file)
#'
#'@export
liftover_junc_id <- function(junc_df, chain_file){


  # convert junc_id into GenomicRanges object
  gr <- junc_to_gr(junc_df$junc_id)

  # read chain file
  chain <- rtracklayer::import.chain(chain_file)

  # Lift over junc_id to genomic position in the personalized genome
  junc_id_lifted_grl  <- rtracklayer::liftOver(gr, chain)

  junc_df %>%
    dplyr::mutate(

      # add lifted junctions ranges as a list and covert back to junc_id
      junc_id_lifted_lst = as.list(junc_id_lifted_grl) %>%
        purrr::map(as.character),

      # check if liftOver was successful and unique
      junc_id_lifted_unique = purrr::map_lgl(junc_id_lifted_lst, ~ length(.x) == 1),

      # collapse multiple IDs with `|` and replace empty junc_ids with NA
      junc_id_lifted = junc_id_lifted_lst %>%
        purrr::map_chr(stringr::str_c, collapse = "|") %>%
        dplyr::na_if("")
    ) %>%

    # remove temporary column
    dplyr::select(-junc_id_lifted_lst)

}




