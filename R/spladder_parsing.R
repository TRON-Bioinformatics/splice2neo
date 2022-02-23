# SplAdder FUNCTION -------------------------------------------------------


#' Sorts columns of junction output file in the following order:
#' "junction_start", "junction_end", "strand", "chromosome", "Gene",
#' "class", "AS_event_ID", "junction_id"
#'
#' @param tib A tibble with the following as given in the description.
#'
#' @return A tibble with sorted columns as given above
#'
#'@import dplyr
#'
sort_columns <- function(tib) {
  column_order <-
    c(
      "junction_start",
      "junction_end",
      "strand",
      "chromosome",
      "Gene",
      "class",
      #"junction_info",
      "AS_event_ID",
      "junc_id"
    )
  tib[column_order]
}


#' Transforms events from alternative 3' splice sites from SplAdder output
#' format into standardized junction format
#'
#' @param tib A tibble in SplAdder output format
#'
#' @return A tibble in standardized junction format
#'
#' @examples
#' spladder_output.a3ss
#' transformed_spladder.a3ss <- spladder_transform_a3ss(spladder_output.a3ss)
#' transformed_spladder.a3ss
#'
#' @import dplyr
#' @export
spladder_transform_a3ss <- function(tib) {
  tib %>%
    mutate(
      junc_id1 = ifelse(
        strand == "+",
        generate_junction_id(contig, exon_const_end, exon_alt1_start, strand),
        generate_junction_id(contig, exon_alt1_end, exon_const_start, strand)
      ),
      junc_id2 = ifelse(
        strand == "+",
        generate_junction_id(contig, exon_const_end, exon_alt2_start, strand),
        generate_junction_id(contig, exon_alt2_end, exon_const_start, strand)
      ),
      junc_id = paste(junc_id1, junc_id2, sep = ";"),
      class = "A3SS",
      AS_event_ID = event_id
    ) %>%
    tidyr::separate_rows(junc_id, sep = ";") %>%
    dplyr::select(junc_id, gene_name, class, AS_event_ID)
}


#' Transforms events from alternative 5' splice sites from SplAdder output
#' format into standardized junction format
#'
#' @param tib A tibble in SplAdder output format
#'
#' @return A tibble in standardized junction format
#'
#' @examples
#'spladder_output.a5ss
#'transformed_spladder.a5ss <- spladder_transform_a5ss(spladder_output.a5ss)
#'transformed_spladder.a5ss
#'
#'@import dplyr
#' @export
spladder_transform_a5ss <- function(tib) {
  tib %>%
    mutate(
      junc_id1 = ifelse(
        strand == "+",
        generate_junction_id(contig, exon_alt1_end, exon_const_start, strand),
        generate_junction_id(contig, exon_const_end, exon_alt1_start, strand)
      ),
      junc_id2 = ifelse(
        strand == "+",
        generate_junction_id(contig, exon_alt2_end, exon_const_start, strand),
        generate_junction_id(contig, exon_const_end, exon_alt2_start, strand)
      ),
      junc_id = paste(junc_id1, junc_id2, sep = ";"),
      class = "A5SS",
      AS_event_ID = event_id
    ) %>%
    tidyr::separate_rows(junc_id, sep = ";") %>%
    dplyr::select(junc_id, gene_name, class, AS_event_ID)
}

#' Transforms events resulting from exon skipping from SplAdder output
#' format into standardized junction format
#'
#' @param tib A tibble in SplAdder output format
#'
#' @return A tibble in standardized junction format
#'
#' @examples
#'spladder_output.exonskip
#'transformed_spladder.exonskip <- spladder_transform_exon_skipping(
#'  spladder_output.exonskip)
#'transformed_spladder.exonskip
#'
#'@import dplyr
#' @export
spladder_transform_exon_skipping <- function(tib) {
  tib %>%
    mutate(
      junc_id1 = generate_junction_id(contig, exon_pre_end, exon_start, strand),
      junc_id2 = generate_junction_id(contig, exon_end, exon_aft_start, strand),
      junc_id3 = generate_junction_id(contig, exon_pre_end, exon_aft_start, strand),
      junc_id = paste(junc_id1, junc_id2, junc_id3, sep = ";"),
      class = "cassette_exon",
      AS_event_ID = event_id
    ) %>%
    separate_rows(junc_id, sep = ";") %>%
    dplyr::select(junc_id, gene_name, class, AS_event_ID)
}


#' Transforms events resulting from intron retention from SplAdder output
#' format into standardized junction format
#'
#' @param tib A tibble in SplAdder output format
#'
#' @return A tibble in standardized junction format
#'
#' @examples
#'spladder_output.intronreten
#'transformed_spladder.intronreten <- spladder_transform_intron_retention(
#'  spladder_output.intronreten)
#'transformed_spladder.intronreten
#'
#'@import dplyr
#' @export
spladder_transform_intron_retention <- function(tib) {
  tib %>%
    mutate(
      junc_id1 = generate_junction_id(contig, exon1_end, intron_start, strand),
      junc_id2 = generate_junction_id(contig, intron_end, exon2_start, strand),
      class = "intron_retention",
      AS_event_ID = event_id,
      junc_id = paste(junc_id1, junc_id2, sep = ";")
    ) %>%
    separate_rows(junc_id, sep = ";") %>%
    dplyr::select(junc_id, gene_name, class, AS_event_ID)
}


#' Transforms events resulting from mutually exclusive exons from SplAdder output
#' format into standardized junction format
#'
#' @param tib A tibble in SplAdder output format
#'
#' @return A tibble in standardized junction format
#'
#' @examples
#'spladder_output.mutexon
#'transformed_spladder.mutexon <- spladder_transform_mutex_exon(
#'  spladder_output.mutexon)
#'transformed_spladder.mutexon
#'
#'@import dplyr
#' @export
spladder_transform_mutex_exon <- function(tib) {
  tib %>%
    mutate(
      junc_id1 = generate_junction_id(contig, exon_pre_end, exon1_start, strand),
      junc_id2 = generate_junction_id(contig, exon1_end, exon_aft_start, strand),
      junc_id3 = generate_junction_id(contig, exon_pre_end, exon2_start, strand),
      junc_id3 = generate_junction_id(contig, exon2_end, exon_aft_start, strand),
      junc_id = paste(junc_id1, junc_id2, junc_id3, sep=";"),
      class = "mutex_exon",
      AS_event_ID = event_id
    ) %>%
    separate_rows(junc_id, sep = ";") %>%
    dplyr::select(junc_id, gene_name, class, AS_event_ID)
}


#' Transforms SplAdder output into standardized junction format
#'
#' @param l A list with tibbles that contain the SplAdder output - each for one type of alternative splicing. These types can be "A5SS",
#'  "A3SS", "cassette_exon", "intron_retention", "mutex_exons".
#'
#' @return A tibble in standardized junction format, combining all alternative
#'   splicing classes that are were determined with SplAdder
#'
#' @examples
#'spladder_output
#'transformed_spladder <- spladder_transform_format(
#'  spladder_output)
#'transformed_spladder
#'
#'@import dplyr
#' @export
spladder_transform_format <- function(l) {

  l_new <- l

  if("A5SS" %in% names(l)){
    l_new$A5SS <- spladder_transform_a5ss(l$A5SS)
  }
  if("A3SS" %in% names(l)){
    l_new$A3SS <- spladder_transform_a3ss(l$A3SS)
  }
  if("cassette_exon" %in% names(l)){
    l_new$cassette_exon <- spladder_transform_exon_skipping(l$cassette_exon)
  }
  if("intron_retention" %in% names(l)){
    l_new$intron_retention <- spladder_transform_intron_retention(l$intron_retention)
  }
  if("mutex_exons" %in% names(l)){
    l_new$mutex_exons <- spladder_transform_mutex_exon(l$mutex_exons)
  }

  df <- bind_rows(
    l_new
  )
  df %>%
    tidyr::separate(
      junc_id,
      into = c("chromosome", "junction_start_end", "strand"),
      sep = ":",
      remove = F
    ) %>%
    tidyr::separate(
      junction_start_end,
      into = c("junction_start", "junction_end"),
      sep = "-",
      remove = F
    ) %>%
    mutate(
      junction_start = as.numeric(junction_start),
      junction_end = as.numeric(junction_end)
    ) %>%
    dplyr::rename(., Gene = gene_name) %>%
    sort_columns()

}


#' Imports SplAdder output from a given path with ".confirmed.txt.gz" files.
#'
#' @param path The path to a folder with SplAdder output. This folder must contain the ".confirmed.txt.gz" files for the alternative splicing type of interest.
#'
#' @return A list with tibbles. Each tibble is a SplAdder output for "A5SS",
#'  "A3SS", "cassette_exon", "intron_retention", "mutex_exons".
#'
#'
#' @import readr
#' @export
import_spladder <- function(path){

  message("Importing Spladder files ...")

  files <- list.files(path, "confirmed.txt.gz")
  as_types <- gsub("_C[0-3].confirmed.txt.gz", "", files)
  as_types <- gsub("merge_graphs_", "", as_types)
  path_files <- paste(path, files ,sep = "/" )

  files <- lapply(path_files, read_delim,
                  delim = "\t",
                  col_types = cols(
                    contig = col_character()
                  ),
                  show_col_types = FALSE)

  if(length(files) == 0){
    stop("There are no SplAdder confirmed.txt.gz input files")
  }
  names_events <-
    c(
      "alt_3prime" = "A3SS",
      "alt_5prime" = "A5SS",
      "exon_skip" = "cassette_exon",
      "intron_retention" = "intron_retention",
      "mutex_exons" = "mutex_exons"
    )
  names(files) <- names_events[as_types]

  message("Importing Spladder files completed")

  return(files)
}


#' Imports SplAdder output from a given path and transforms it into standardized
#' junction format
#'
#' @param path The path to a folder with SplAdder output
#'
#' @return A tibble in standardized junction format, combining all alternative
#'   splicing classes that are covered by SplAdder
#'
#'
#' @import readr
#' @export
spladder_transform <- function(path){
  dat <- import_spladder(path)
  juncs <- spladder_transform_format(dat)
  return(juncs)
}
