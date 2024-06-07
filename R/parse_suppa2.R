# SUPPA2 FUNCTIONS -----------------------------------------------------

#' Imports "_strict.ioe" from SUPPA2 output.
#'
#' @param file.name The path to the file
#'
#' @return A tibble with columns "seqnames", "gene_id" "event_id"
#'         "inclusion_transcripts" and "total_transcripts"
#'
#'
#' @import readr
#' @export
read_suppa_ioe <- function(file.name) {
  if(!file.exists(file.name)){
    stop("SUPPA ioe files are missing")
  }
  dat_suppa <- read_delim(file.name,
    delim = "\t",
    col_names = T,
    show_col_types = FALSE,
    col_types = cols(.default = "c")
  )
  return(dat_suppa)
}


#' Transforms events resulting from exon skipping from SUPPA2 output
#' format into standardized junction format
#'
#' @param tib A tibble in SUPPA2 output format
#'
#' @return A tibble in standardized junction format
#'
#'
#'@import dplyr tidyr
transform_suppa_se_events <- function(tib) {
  tib <- tib %>%
    tidyr::separate(
      event_id,
      sep = ":",
      into = c("gene", "chromosome", "junc1", "junc2", "strand"),
      remove = F
    ) %>%
    tidyr::separate(
      junc1,
      sep = "-",
      into = c("e1_end", "e2_start")
    ) %>%
    tidyr::separate(
      junc2,
      sep = "-",
      into = c("e2_end", "e3_start")
    ) %>%
    dplyr::mutate(
      Gene = str_replace(gene,";SE", ""),
      class = "cassette_exon",
      AS_event_ID = event_id,
      junc_id1 = generate_junction_id(chromosome, e1_end, e2_start, strand),
      junc_id2 = generate_junction_id(chromosome, e2_end, e3_start, strand),
      junc_id3 = generate_junction_id(chromosome, e1_end, e3_start, strand),
      junc_id = paste(junc_id1, junc_id2, junc_id3, sep = ";")
    ) %>%
    tidyr::separate_rows(junc_id, sep=";") %>%
    dplyr::select(junc_id, Gene, class, AS_event_ID)
  return(tib)
}

#' Transforms events resulting from intron retention from SUPPA2 output
#' format into standardized junction format
#'
#' @param tib A tibble in SUPPA2 output format
#'
#' @return A tibble in standardized junction format
#'
#'
#'@import dplyr tidyr
transform_suppa_ir_events <- function(tib) {
  tib <- tib %>%
    tidyr::separate(
      event_id,
      sep = ":",
      into = c("gene", "chromosome", "e1_start", "intron", "e2_end", "strand"),
      convert = TRUE,
      remove = F
    ) %>%
    tidyr::separate(
      intron,
      sep = "-",
      into = c("e1_end", "e2_start"),
      convert = TRUE
    ) %>%
    dplyr::mutate(
      Gene = str_replace(gene,";RI", ""),
      class = "intron_retention",
      AS_event_ID = event_id,
      junc_id1 = generate_junction_id(chromosome, e1_end, e1_end + 1, strand),
      junc_id2 = generate_junction_id(chromosome, e2_start - 1 , e2_start, strand),
      junc_id = paste(junc_id1, junc_id2, sep = ";")
    ) %>%
    tidyr::separate_rows(junc_id, sep=";") %>%
    dplyr::select(junc_id, Gene, class, AS_event_ID)
  return(tib)
}

#' Transforms events from alternative 3' or 5' splice sites from SUPPA2 output
#' format into standardized junction format
#'
#' @param tib A tibble in SUPPA2 output format
#' @param type The type of alternative splice site. Can be "A3SS" or "A5SS"
#'
#' @return A tibble in standardized junction format
#'
#'
#' @import dplyr tidyr
transform_suppa_ass_events <- function(tib, type="A5SS") {
  stopifnot(
    type %in% c("A5SS","A3SS"),
    "{type} not of type A3SS or A5SS")
  tib <- tib %>%
    tidyr::separate(
      event_id,
      sep = ":",
      into = c("gene", "chromosome", "junc1", "junc2", "strand"),
      remove = F
    ) %>%
    # junc1 is either e2_end and e3_start for A5SS:+ and A3SS:-
    # or e1_end and e2_start for A5SS:- and A3SS:+
    tidyr::separate(
      junc1,
      sep = "-",
      into = c("junc1_start", "junc1_end")
    ) %>%
    # junc2 always consists of e1_end & e3_start
    # no matter which event type or strand were are on
    tidyr::separate(
      junc2,
      sep = "-",
      into = c("e1_end", "e3_start")
    ) %>%
    dplyr::mutate(
      Gene = case_when(
        type == "A5SS" ~ str_replace(gene, ";A5", ""),
        type == "A3SS" ~ str_replace(gene, ";A3", "")
      ),
      class = type,
      AS_event_ID = event_id,
      junc_id1 = case_when(
        type == "A5SS" & strand == "+" | type == "A3SS" & strand == "-" ~ generate_junction_id(chromosome, junc1_start, e3_start, strand),
        type == "A5SS" & strand == "-" | type == "A3SS" & strand == "+" ~ generate_junction_id(chromosome, e1_end, junc1_end, strand)
      ),
      junc_id2 = generate_junction_id(chromosome, e1_end, e3_start, strand),
      junc_id = paste(junc_id1, junc_id2, sep = ";")
    ) %>%
    tidyr::separate_rows(junc_id, sep=";") %>%
    dplyr::select(junc_id, Gene, class, AS_event_ID)

  return(tib)
}

#' Transforms events resulting from mutually exclusive exons from SUPPA2 output
#' format into standardized junction format
#'
#' @param tib A tibble in SUPPA2 output format
#'
#' @return A tibble in standardized junction format
#'
#'
#' @import dplyr tidyr
transform_suppa_mxe_events <- function(tib) {
  tib <- tib %>%
    tidyr::separate(
      event_id,
      sep = ":",
      into = c("gene", "chromosome", "junc1", "junc2", "junc3", "junc4", "strand"),
      remove = F
    ) %>%
    dplyr::mutate(
      e1_end = str_split_fixed(junc1, "-", n = 2)[, 1],
      e2_start = str_split_fixed(junc1, "-", n = 2)[, 2],
      e2_end = str_split_fixed(junc2, "-", n = 2)[, 1],
      e3_start = str_split_fixed(junc3, "-", n = 2)[, 2],
      e3_end = str_split_fixed(junc4, "-", n = 2)[, 1],
      e4_start = str_split_fixed(junc4, "-", n = 2)[, 2],
    ) %>%
    dplyr::mutate(
      Gene = str_replace(gene,";MX", ""),
      class = "mutex_exons",
      AS_event_ID = event_id,
      junc_id1 = generate_junction_id(chromosome, e1_end, e2_start, strand),
      junc_id2 = generate_junction_id(chromosome, e2_end, e4_start, strand),
      junc_id3 = generate_junction_id(chromosome, e1_end, e3_start, strand),
      junc_id4 = generate_junction_id(chromosome, e3_end, e4_start, strand),
      junc_id = paste(junc_id1, junc_id2, junc_id3, junc_id4, sep = ";")
    ) %>%
    tidyr::separate_rows(junc_id, sep=";") %>%
    dplyr::select(junc_id, Gene, class, AS_event_ID)
  return(tib)
}

#' Transforms SUPPA2 ioe event files file into standardized junction format.
#'
#' @param tib_list A list with tibbles that contain the SUPPA2 output - each for one type of alternative splicing.
#'    These types can be "A5SS", "A3SS", "cassette_exon", "intron_retention", "mutex_exons".
#'
#' @return A tibble in standardized junction format, combining all alternative
#'   splicing classes that are were determined with SUPPA2
#'
#'
#' @import dplyr
#' @export
suppa_transform_format <- function(tib_list) {
  l_new <- tib_list

  if("A5SS" %in% names(tib_list)){
    l_new$A5SS <- transform_suppa_ass_events(l_new$A5SS, type = "A5SS")
  }
  if("A3SS" %in% names(tib_list)){
    l_new$A3SS <- transform_suppa_ass_events(l_new$A3SS, type = "A3SS")
  }
  if("cassette_exon" %in% names(tib_list)){
    l_new$cassette_exon <- transform_suppa_se_events(l_new$cassette_exon)
  }
  if("intron_retention" %in% names(tib_list)){
    l_new$intron_retention <- transform_suppa_ir_events(l_new$intron_retention)
  }
  if("mutex_exons" %in% names(tib_list)){
    l_new$mutex_exons <- transform_suppa_mxe_events(l_new$mutex_exons)
  }

  df <- dplyr::bind_rows(
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
    dplyr::select(
      chromosome,
      junction_start,
      junction_end,
      junc_id,
      strand,
      Gene,
      class,
      AS_event_ID
    ) %>%
    sort_columns() %>%
    dplyr::distinct(junc_id, .keep_all = TRUE)

}

#' Imports SUPPA2 output from a given path with "_strict.ioe" files.
#' The results for one patient should be stored in the given path. Please note
#' that AF and AL events are not supported.
#'
#' @param path The path to a folder with SUPPA2 output. This folder must contain the "_strict.ioe"
#' files for the alternative splicing type of interest.
#'
#' @return A list with tibbles. Each tibble is a SUPPA2 output for "A5SS",
#'  "A3SS", "cassette_exon", "intron_retention", "mutex_exons".
#'
#'
#' @import readr purrr
#' @export
suppa_import <- function(path) {
  message("Importing SUPPA files ...")

  files <- list.files(path, "_strict.ioe")
  if(any(grepl("AL_strict.ioe", files))){
    files <- files[-which(grepl("AL_strict.ioe", files))]
  }
  if(any(grepl("AF_strict.ioe", files))){
    files <- files[-which(grepl("AL_strict.ioe", files))]
  }
  as_types <- gsub("_strict.ioe", "", files)
  n_items <- length(str_split(as_types, "_")[[1]])
  as_types <- sapply(str_split(as_types, "_"), "[[", n_items)
  path_files <- paste(path, files ,sep = "/" )

  files <- purrr::map(path_files, read_suppa_ioe)

  if(length(files) == 0){
    stop("There are no SUPPA ioe input files")
  }
  names_events <-
    c(
      "A3" = "A3SS",
      "A5" = "A5SS",
      "SE" = "cassette_exon",
      "RI" = "intron_retention",
      "MX" = "mutex_exons"
    )
  names(files) <- names_events[as_types]
  message("Importing SUPPA files completed")
  return(files)
}

#' Imports SUPPA2 output from a given path and transforms it into standardized
#' junction format. Results from one patients should be stored per folder. Please note
#' that alternative first exons (AF) and alternative last exons (AL) are not supported.
#'
#' - GitHub: https://github.com/comprna/SUPPA
#' - Paper: https://doi.org/10.1186/s13059-018-1417-1
#'
#' @param path The path to a folder with SUPPA2 output
#'
#' @return A tibble in standardized junction format, combining all alternative
#'   splicing classes that are covered by SUPPA2
#'
#' @import readr
#' @export
suppa_transform <- function(path) {
  dat <- suppa_import(path)
  juncs <- suppa_transform_format(dat)
  return (juncs)
}
