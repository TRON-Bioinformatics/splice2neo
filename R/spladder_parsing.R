# SplAdder FUNCTION -------------------------------------------------------


#' Sorts columns of junction output file in the following order:
#' "junction_start", "junction_end", "strand", "chromosome", "Gene",
#' "class", "junction_id"
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
      "junc_id",
      "number_supporting_reads"
    )
  tib[column_order]
}


#' Transforms events from alternative 3' or 5' splice sites from SplAdder output
#' format into standardized junction format
#'
#' @param tib A tibble in SplAdder output format
#' @param AS_class The type of alternative splice site. Can be "A3SS" or "A5SS"
#'
#' @return A tibble in standardized junction format
#'
#'
#' @import dplyr
spladder_transform_ass <- function(tib, AS_class) {
  tib %>%
    mutate(
      junc_id1 = ifelse(
        strand == "+" & AS_class == "A3SS" | strand == "-" & AS_class == "A5SS",
        generate_junction_id(chrm, e1_end, e2_start, strand),
        generate_junction_id(chrm, e1_end, e3_start, strand)
      ),
      junc_id2 = ifelse(
        strand == "+" & AS_class == "A3SS" | strand == "-" & AS_class == "A5SS",
        generate_junction_id(chrm, e1_end, e3_start, strand),
        generate_junction_id(chrm, e2_end, e3_start, strand)
      ),
      supporting_reads_junc1 = ifelse(
        strand == "+" & AS_class == "A3SS" | strand == "-" & AS_class == "A5SS",
        get(grep('e2_conf', colnames(tib), value = T)),
        get(grep('e1e3_conf', colnames(tib), value = T))
      ),
      supporting_reads_junc2 = ifelse(
        strand == "+" & AS_class == "A3SS" | strand == "-" & AS_class == "A5SS",
        get(grep('e1e3_conf', colnames(tib), value = T)),
        get(grep('e2_conf', colnames(tib), value = T))
      ),
      junc_id = paste(junc_id1, junc_id2, sep = ";"),
      number_supporting_reads = paste(supporting_reads_junc1, supporting_reads_junc2, sep = ";"),
      class = AS_class,
    ) %>%
    tidyr::separate_rows(junc_id, number_supporting_reads, sep = ";") %>%
    dplyr::select(junc_id, gene_name, class, AS_event_ID, number_supporting_reads) %>%
    dplyr::mutate(number_supporting_reads = as.numeric(number_supporting_reads))
}



#' Transforms events resulting from exon skipping from SplAdder output
#' format into standardized junction format
#'
#' @param tib A tibble in SplAdder output format
#'
#' @return A tibble in standardized junction format
#'
#'
#'@import dplyr
spladder_transform_exon_skipping <- function(tib) {
  tib %>%
    mutate(
      junc_id1 = generate_junction_id(chrm, e1_end, e2_start, strand),
      junc_id2 = generate_junction_id(chrm, e2_end, e3_start, strand),
      junc_id3 = generate_junction_id(chrm, e1_end, e3_start, strand),
      supporting_reads_junc1 = get(grep('e1e2_conf', colnames(tib), value = T)),
      supporting_reads_junc2 = get(grep('e2e3_conf', colnames(tib), value = T)),
      supporting_reads_junc3 = get(grep('e1e3_conf', colnames(tib), value = T)),
      junc_id = paste(junc_id1, junc_id2, junc_id3, sep = ";"),
      number_supporting_reads = paste(supporting_reads_junc1,
                                     supporting_reads_junc2,
                                     supporting_reads_junc3, sep = ";"),
      class = "cassette_exon",
      AS_event_ID = event_id
    ) %>%
    tidyr::separate_rows(junc_id, number_supporting_reads, sep = ";") %>%
    dplyr::select(junc_id, gene_name, class, AS_event_ID, number_supporting_reads) %>%
    dplyr::mutate(number_supporting_reads = as.numeric(number_supporting_reads))
}


#' Transforms events resulting from intron retention from SplAdder output
#' format into standardized junction format
#'
#' @param tib A tibble in SplAdder output format
#'
#' @return A tibble in standardized junction format
#'
#'
#'@import dplyr
spladder_transform_intron_retention <- function(tib) {
  tib %>%
    mutate(
      junc_id1 = generate_junction_id(chrm, e1_end, e2_start, strand),
      junc_id2 = generate_junction_id(chrm, e2_end, e3_start, strand),
      number_supporting_reads = get(grep('e1e3_conf', colnames(tib), value = T)),
      class = "intron_retention",
      junc_id = paste(junc_id1, junc_id2, sep = ";")
    ) %>%
    tidyr::separate_rows(junc_id, sep = ";") %>%
    dplyr::select(junc_id, gene_name, class, AS_event_ID, number_supporting_reads) %>%
    dplyr::mutate(number_supporting_reads = as.numeric(number_supporting_reads))
}


#' Transforms events resulting from mutually exclusive exons from SplAdder output
#' format into standardized junction format
#'
#' @param tib A tibble in SplAdder output format
#'
#' @return A tibble in standardized junction format
#'
#'
#'@import dplyr
spladder_transform_mutex_exon <- function(tib) {
  tib %>%
    mutate(
      junc_id1 = generate_junction_id(chrm, e1_end, e2_start, strand),
      junc_id2 = generate_junction_id(chrm, e2_end, e4_start, strand),
      junc_id3 = generate_junction_id(chrm, e1_end, e3_start, strand),
      junc_id4 = generate_junction_id(chrm, e3_end, e4_start, strand),
      supporting_reads_junc1 = get(grep('e1e2_conf', colnames(tib), value = T)),
      supporting_reads_junc2 = get(grep('e2e4_conf', colnames(tib), value = T)),
      supporting_reads_junc3 = get(grep('e1e3_conf', colnames(tib), value = T)),
      supporting_reads_junc4 = get(grep('e3e4_conf', colnames(tib), value = T)),
      junc_id = paste(junc_id1, junc_id2, junc_id3, junc_id4, sep=";"),
      number_supporting_reads = paste(supporting_reads_junc1,
                                     supporting_reads_junc2,
                                     supporting_reads_junc3,
                                     supporting_reads_junc4, sep = ";"),
      class = "mutex_exon",
      AS_event_ID = event_id
    ) %>%
    tidyr::separate_rows(junc_id, number_supporting_reads, sep = ";") %>%
    dplyr::select(junc_id, gene_name, class, AS_event_ID, number_supporting_reads) %>%
    dplyr::mutate(number_supporting_reads = as.numeric(number_supporting_reads))
}


#' Transforms SplAdder output into standardized junction format
#'
#' @param l A list with tibbles that contain the SplAdder output - each for one type of alternative splicing. These types can be "A5SS",
#'  "A3SS", "cassette_exon", "intron_retention", "mutex_exons".
#'
#' @return A tibble in standardized junction format, combining all alternative
#'   splicing classes that are were determined with SplAdder
#'
#'@import dplyr
#' @export
spladder_transform_format <- function(l) {

  l_new <- l

  if("A5SS" %in% names(l)){
    l_new$A5SS <- spladder_transform_ass(l$A5SS, AS_class = "A5SS")
  }
  if("A3SS" %in% names(l)){
    l_new$A3SS <- spladder_transform_ass(l$A3SS, AS_class = "A3SS")
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
    )%>%
    dplyr::rename(., Gene = gene_name) %>%
    sort_columns()%>%
    distinct(junc_id, .keep_all = TRUE)

}


#' Imports SplAdder output from a given path with ".confirmed.txt.gz" files.
#' The results for one patient should be stored in the given path. Please note
#' that multiple (coordinated) exons skips (mult_exon_skip) are currently not supported.
#'
#' @param path The path to a folder with SplAdder output. This folder must contain the ".confirmed.txt.gz"
#' files for the alternative splicing type of interest.
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
  if(any(grepl("mult_exon_skip", files))){
    files <- files[-which(grepl("mult_exon_skip", files))]
  }
  as_types <- gsub("_C[0-3].confirmed.txt.gz", "", files)
  str_split(as_types, "_")
  n_items <- length(str_split(as_types, "_")[[1]])
  typ1 <- sapply(str_split(as_types, "_"), "[[", n_items-1)
  typ2 <- sapply(str_split(as_types, "_"), "[[", n_items)
  as_types <- paste0(typ1, "_", typ2)
  path_files <- paste(path, files ,sep = "/" )

  files <- lapply(path_files, read_delim,
                  delim = "\t",
                  col_types = cols(.default = "c"))

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
#' junction format. Results from one patients should be stored per folder.
#'
#' @param path The path to a folder with SplAdder output
#'
#' @return A tibble in standardized junction format, combining all alternative
#'   splicing classes that are covered by SplAdder
#' @examples
#' path <-  system.file("extdata", "", package = "splice2neo")
#' spladder_juncs <- spladder_transform(path)
#' spladder_juncs
#'
#' @import readr
#' @export
spladder_transform <- function(path){
  dat <- import_spladder(path)
  juncs <- spladder_transform_format(dat)
  return(juncs)
}
