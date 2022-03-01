
# COMBINED DATASET --------------------------------------------------------


#' Combines tibbles with junctions from SplAdder and LeafCutter into a combined dataset of expressed splice junctions
#'
#' @param leafcutter_juncs A tibble the junctions identified by LeafCutter in
#'   standardized format
#' @param spladder_juncs A tibble the junctions identified by SplAdder in
#'   standardized format
#'
#' @return A combined table with unique junctions. The columns RNA_tool
#'   contains information which tools identified the given junction
#'
#' @examples
#' path <-  system.file("extdata", "", package = "splice2neo")
#' spladder_juncs <- spladder_transform(path)
#' leafcutter_juncs <- leafcutter_transform(path)
#' dat.combined <- generate_combined_dataset(spladder_juncs, leafcutter_juncs)
#'
#' @import dplyr
#' @export
generate_combined_dataset <- function(spladder_juncs, leafcutter_juncs){

  rna_juncs <- spladder_juncs %>%
    bind_rows(leafcutter_juncs) %>%
    distinct(junc_id, .keep_all = T) %>%
    mutate(
      identified_by_leafcutter = ifelse(junc_id %in% leafcutter_juncs$junc_id, TRUE, FALSE)
      ,
      identified_by_spladder = ifelse(junc_id %in% spladder_juncs$junc_id, TRUE, FALSE)
    )
  return(rna_juncs)

}

#' This is a wrapper function to directly map the information if a junction predicted from WES data was found in RNA-seq by LeafCutter or SplAdder
#'
#' @param path_to_spladder The path to the results from RNA-seq analysis with SplAdder
#' @param path_to_leafcutter The path to the results from RNA-seq analysis with LeafCutter
#' @param mutation_juncs The junction-transcript centric data predicted from WES data.
#'
#' @return The junction-transcript centric data predicted from WES data is extended by the information if a respective aberrant junctions was
#' identified by spladder or LeafCutter (`identified_by_leafcutter`, `identified_by_spladder`)
#'
#' @examples
#'
#' @import dplyr
#' @export
#'
add_identified_in_RNA <- function(mutation_juncs, path_to_spladder, path_to_leafcutter){

  spladder_juncs <- spladder_transform(path_to_spladder)
  leafcutter_juncs <- leafcutter_transform(path_to_leafcutter)

  rna_juncs <- generate_combined_dataset(spladder_juncs = spladder_juncs, leafcutter_juncs = leafcutter_juncs)

  rna_juncs <- rna_juncs %>%
    select(junc_id, identified_by_leafcutter, identified_by_spladder)

  mutation_juncs <- mutation_juncs %>%
    left_join(rna_juncs, by = "junc_id")

  return(mutation_juncs)

}
