
# COMBINED DATASET --------------------------------------------------------


#' Combines tibbles with junctions from any number of RNA-seq tools into a
#' combined dataset of expressed splice junctions
#'
#' @param rna_junc_data_list A named list with junction tibbles in standardized
#'   format.
#'
#' @return A combined table with unique junctions. The columns
#'   identified_by_{name} contains information which tools identified the given
#'   junction
#'
#' @examples
#' path <-  system.file("extdata", "", package = "splice2neo")
#' spladder_juncs <- spladder_transform(path)
#' path <-  system.file("extdata", "test_regtools_Aligned.out.sorted.bam.junc", package = "splice2neo")
#' regtools_juncs <- regtools_transform(path)
#' dat.combined <- generate_combined_dataset(list("spladder" = spladder_juncs,
#'                                                "regtools" = regtools_juncs))
#'
#' @import dplyr purrr stringr tidyr
#' @export
generate_combined_dataset <- function(rna_junc_data_list){

  stopifnot(is.list(rna_junc_data_list))
  stopifnot(!is.null(names(rna_junc_data_list)))
  stopifnot(all(map_lgl(rna_junc_data_list, is.data.frame)))
  stopifnot(all(map_lgl(rna_junc_data_list, ~'junc_id' %in% colnames(.x))))
  stopifnot(all(map_lgl(rna_junc_data_list, ~'junction_start' %in% colnames(.x))))
  stopifnot(all(map_lgl(rna_junc_data_list, ~'junction_end' %in% colnames(.x))))
  stopifnot(all(map_lgl(rna_junc_data_list, ~'strand' %in% colnames(.x))))
  stopifnot(all(map_lgl(rna_junc_data_list, ~'chromosome' %in% colnames(.x))))

  indicator_columns <-
    c("junc_id", "junction_start", "junction_end", "strand", "chromosome")

  rna_juncs <- rna_junc_data_list %>%
    dplyr::bind_rows(.id = "tool") %>%
    dplyr::mutate(detected = TRUE) %>%
    dplyr::select(all_of(c(indicator_columns, "tool", "detected"))) %>%
    tidyr::complete(
      nesting(junc_id, junction_start, junction_end, strand, chromosome),
        tool, fill = list(detected = FALSE)) %>%
    tidyr::pivot_wider(names_from = "tool", values_from = "detected",
        names_glue = "identified_by_{.name}")

  # rename annotation columns
  my_rename <- function(x, prefix_name){
    stringr::str_c(prefix_name, "_", x)
  }

  # for each input data add the tool/source name to all columm names
  # except the indicator columns
  rna_junc_data_list_names <-
      purrr::map2(rna_junc_data_list,
                  names(rna_junc_data_list),
                  ~rename_with(
                      .data = .x,
                      .fn = my_rename,
                      .cols = -all_of(indicator_columns),
                      prefix_name = .y))

  # add tool specific columns/annotations
  # by iteratively apply a left_join()
  for (this_df in rna_junc_data_list_names){
    rna_juncs <- rna_juncs %>%
      dplyr::left_join(this_df,
          by = indicator_columns)
  }
  # Rename shared columns that are not indicator columns. For now they are collapsed by comma
  rna_juncs <- rna_juncs %>%
    tidyr::unite("Gene", ends_with("_Gene"), sep = ",", remove = TRUE, na.rm = TRUE)
  rna_juncs <- rna_juncs %>%
    tidyr::unite("class", ends_with("_class"), sep=",", remove = TRUE, na.rm = TRUE)
  rna_juncs <- rna_juncs %>%
    tidyr::unite("AS_event_ID", ends_with("_AS_event_ID"), sep = ",", remove = TRUE, na.rm = TRUE)

  return(rna_juncs)
}



#' This is a wrapper function to directly map the information if a junction predicted from WES data was found in RNA-seq by Regtools or SplAdder
#'
#' @param path_to_spladder The path to the results from RNA-seq analysis with SplAdder
#' @param path_to_regtools The path to the results from RNA-seq analysis with Regtools
#' @param mutation_juncs The junction-transcript centric data predicted from WES data.
#'
#' @return The junction-transcript centric data predicted from WES data is extended by the information if a respective aberrant junctions was
#' identified by spladder or regtools (`identified_by_regtools`, `identified_by_spladder`)
#'
#'
#' @import dplyr
#' @export
#'
add_identified_in_RNA <- function(mutation_juncs, path_to_spladder, path_to_regtools){

  spladder_juncs <- spladder_transform(path_to_spladder)
  regtools_juncs <- regtools_transform(path_to_regtools)

  rna_juncs <- generate_combined_dataset(list("spladder" = spladder_juncs, "regtools" = regtools_juncs))

  rna_juncs <- rna_juncs %>%
    select(junc_id, identified_by_regtools, identified_by_spladder)

  mutation_juncs <- mutation_juncs %>%
    left_join(rna_juncs, by = "junc_id")

  return(mutation_juncs)

}
