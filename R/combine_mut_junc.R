#' Combines data sets with junctions from several different sources
#'
#' @param junc_data_list A named list of tibbles with junctions. the name
#'    should be the source, e.g. the tool name, such as `spliceai` or `pangolin`.
#'    The individual tibbles should at least contain the columns `mut_id`,
#'    `tx_id`, `junc_id` and might have individual sets of other tool/source
#'    specific columns.
#'
#' @return A combined data set with unique junctions based on the four columns
#'    `mut_id`, `tx_id`, `junc_id`, `event_type`.
#'    Additional columns in the input data.frames will be prefixed with the
#'    tool/source name followed by an underscore `_`.
#'    E.g. the column `score` in the input data sets `spliceai` becomes
#'    `spliceai_score`. Additionally, for each tool/source with name `<name>`
#'    an additional column `<name>_detected` will be added to indicate if a
#'    given junction was detected with the indicated tool.
#'
#' @import dplyr
#' @export
combine_mut_junc <- function(junc_data_list){

  # check input types and format
  stopifnot(is.list(junc_data_list))
  stopifnot(!is.null(names(junc_data_list)))
  stopifnot(all(purrr::map_lgl(junc_data_list, is.data.frame)))
  stopifnot(all(purrr::map_lgl(junc_data_list, ~ "mut_id" %in% names(.x))))
  stopifnot(all(purrr::map_lgl(junc_data_list, ~ "junc_id" %in% names(.x))))
  stopifnot(all(purrr::map_lgl(junc_data_list, ~ "tx_id" %in% names(.x))))

  # get union of junctions from all inputs with detection variables
  junc_df <- junc_data_list %>%

    # select distinct junctions by the indicator columns
    purrr::map(dplyr::distinct, mut_id, tx_id, junc_id, event_type) %>%

    # combine into a single data.frame with tool column
    dplyr::bind_rows(.id = "tool") %>%

    # mark as detected
    mutate(detected = TRUE) %>%

    # expand by junction and tool
    complete(
      nesting(mut_id, tx_id, junc_id, event_type),
      tool,
      fill = list(detected = FALSE)
    ) %>%

    # add tools as separate columns
    pivot_wider(names_from = "tool", values_from = "detected",
                names_glue = "{.name}_detected")


  # rename annotation columns
  my_rename = function(x, prefix_name){
    stringr::str_c(prefix_name, "_", x)
  }

  # for each input data add the suffix of the tool/source name to all column names
  # except the indicator columns
  junc_data_list_names = purrr::map2(junc_data_list, names(junc_data_list),
                              ~rename_with(
                                .data = .x,
                                .fn = my_rename,
                                .cols = -all_of(c("mut_id", "tx_id", "junc_id", "event_type")),
                                prefix_name = .y
                              ))

  # add tool/source specific columns/annotations
  # by iteratively apply a left_join()

  for (df in junc_data_list_names){
    junc_df <- junc_df %>%
      left_join(df, by = c("mut_id", "tx_id", "junc_id", "event_type"))
  }

  return(junc_df)

}

