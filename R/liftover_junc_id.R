

#' LiftOver junction IDs using the liftOver tool
#'
#' This can be used to liftOver junctions to personalized genome coordinates.
#'
#' @param junc_df a data.frame with at least one column `junc_id` containing junction IDs
#' @param chain_file path to a chain file for the UCSC liftOver tool. See also \code{\link[rtracklayer]{liftOver}}
#' @return a data.frame like the input `junc_df` with the following additional columns:
#'   - `liftover_successful` a logical vector indicating if the liftOver was successful and lifted junction positions build valid genomic intervals and not single positions.
#'   - `liftover_unique` a logical vector indicating if the liftOver was unique (1-to-1 correspondence).
#'   - `junc_id_lifted_collapsed` a character vector with the lifted junction IDs.
#'   Multiple IDs are separated by `|`.
#'   NA represent junc_ids that could not be lifted.
#'   - `junc_id_lifted` a character vector with a unique lifted junction IDs.
#'   Potentially multiple lifted IDs are combined by the minimal start and maximal
#'   end coordinate. NA represent junc_ids that could not be lifted.
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

  # Lift over junc ranges to genomic position in the personalized genome
  junc_id_lifted_grl  <- rtracklayer::liftOver(gr, chain)

  junc_df %>%
    dplyr::mutate(

      # add lifted junctions ranges as a list
      junc_id_lifted_lst = as.list(junc_id_lifted_grl),

      # check that lifted junctions are valid intervals (not single positions)
      liftover_successful = purrr::map_lgl(junc_id_lifted_lst, ~ length(.x) > 0 & any(BiocGenerics::width(.x) >= 2)),

      # check if liftOver was unique
      liftover_unique = purrr::map_lgl(junc_id_lifted_lst, ~ length(.x) == 1),

      # covert back to junc_id, collapse multiple IDs with `|`, and replace empty junc_id with NA
      junc_id_lifted_collapsed = junc_id_lifted_lst %>%
        purrr::map(as.character) %>%
        purrr::map_chr(stringr::str_c, collapse = "|") %>%
        dplyr::na_if(""),

      # build unique lifted junctions by minimal start and maximal end
      # coordinate of potentially multiple lifted junction ranges
      junc_id_lifted = junc_id_lifted_lst %>%

        # take min start and max end across potentially multiple junctions
        purrr::map(base::range) %>%

        # convert to junc_id and replace empty entries with NA
        purrr::map(as.character) %>%
        purrr::map_if(is_empty, ~ NA_character_) %>%
        unlist()

    ) %>%

    # remove temporary column
    dplyr::select(-junc_id_lifted_lst)

}




