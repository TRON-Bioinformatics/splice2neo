#' Filter junctions for unique combinations of mutation, transcript, and
#' junctions.
#'
#' Filters for unqiue `mut_id`,`junc_id`, `tx_id`. If there are multiple
#' effects (DL/AL) predicted to lead to the same combination of `mut_id`,
#' `junc_id`, `tx_id`, only the  effect with the maximal score is returned.
#'
#' @param junc_data A tibble with junctions. The tibble should contain the columns `mut_id`,
#'    `tx_id`, `junc_id` and `score` and might have individual sets of other tool/source
#'    specific columns. This tibble can be created with `annotate_mut_effect()`.
#'
#' @return A filtered tibble with unique junctions based on the columns
#'    `mut_id`, `tx_id`, `junc_id`. The event with the maximal
#'    `score` is kept.
#'
#' @import dplyr
#' @export
unique_mut_junc <- function(junc_data) {
  # check input types and format
  stopifnot(!is.null(names(junc_data)))
  stopifnot("mut_id" %in% names(junc_data))
  stopifnot("junc_id" %in% names(junc_data))
  stopifnot("tx_id" %in% names(junc_data))
  stopifnot("score" %in% names(junc_data))

  junc_data_res <- junc_data %>%
    arrange(mut_id, junc_id, tx_id, event_type, desc(score)) %>%
    distinct(mut_id, junc_id, tx_id, .keep_all = TRUE)

  return(junc_data_res)

}

#' Filter junction datasets for unique `mut_id`,
#' `junc_id`, `tx_id`. If there were multiple exon inclusion events
#' predicted to lead to the same combination of `mut_id`, `junc_id`, `tx_id`,
#' only the  effect with the maximal score is returned.
#'
#' @param junc_data A tibble with junctions. The tibble should contain the columns `mut_id`,
#'    `tx_id`, `junc_id` and `delta_logit_psi` and might have individual sets of other tool/source
#'    specific columns. This tibble can be created with `annotate_mmsplice()`.
#'
#' @return A filtered tibble with unique junctions based on the columns
#'    `mut_id`, `tx_id`, `junc_id`. The event with the maximal
#'    `delta_logit_psi` is kept.
#'
#' @import dplyr
#' @export
unique_junc_mmsplice <- function(junc_data) {
  # check input types and format
  stopifnot(!is.null(names(junc_data)))
  stopifnot("mut_id" %in% names(junc_data))
  stopifnot("junc_id" %in% names(junc_data))
  stopifnot("tx_id" %in% names(junc_data))
  stopifnot("delta_logit_psi" %in% names(junc_data))

  junc_data_skp <- junc_data %>%
    filter(event_type == "exon skipping")

  junc_data_incl <- junc_data %>%
    filter(event_type == "exon inclusion") %>%
    arrange(mut_id, junc_id, tx_id, desc(delta_logit_psi)) %>%
    distinct(mut_id, junc_id, tx_id, .keep_all = TRUE)

  junc_data_res <- bind_rows(junc_data_skp, junc_data_incl) %>%
    arrange(mut_id, junc_id, tx_id)

  return(junc_data_res)

}
