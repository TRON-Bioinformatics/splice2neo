# LEAFCUTTER FUNCTIONS -----------------------------------------------------

#' Imports "_perind.counts.gz" from LeafCutter output.
#'
#' @param file.name The path to the file
#'
#' @return A tibble with columns "intron_cluster" and "counts
#'
#'
#'
#' @import readr
#' @export
import_leafcutter_counts <- function(file.name) {
  if(!file.exists(file.name)){
    stop("Aligned.out.bam.junc file is missing")
  }
  myData <- read_delim(file.name,
                       delim = " ",
                       skip = 1,
                       col_names = F,
                       show_col_types = FALSE)
  col.nams <- c("intron_cluster", "counts")
  colnames(myData) <- col.nams
  return(myData)
}


#' Transforms LeafCutter counts file into standardized junction format.
#'
#' @param tib LeafCutter counts file as tibble
#'
#' @return A tibble with columns "chromosome", "junction_start", "junction_end",
#'  "cluster"
#'
#'
#' @import dplyr
#' @export
transform_leafcutter_counts <- function(tib) {
  tib %>%
    separate(
      intron_cluster,
      sep = ":",
      into = c("chromosome", "junction_start", "junction_end", "cluster")
    ) %>%
    separate(
      cluster,
      sep = "_",
      into = c("clu", "clu_no", "strand")
    ) %>%
    mutate(
      Gene = NA,
      class = NA,
      AS_event_ID = NA,
      junc_id = generate_junction_id(chromosome, junction_start, junction_end, strand)
    )  %>%
    dplyr::select(
      junction_start,
      junction_end,
      strand,
      chromosome,
      Gene,
      class,
      AS_event_ID,
      junc_id
    )
}






#' Imports "*perind.counts.gz" from LeafCutter output and transforms the raw output
#' into standardized junction output format. Results of one patient should be stored in the given path.
#'
#' @param path The path to leafcutter output
#'
#' @return A tibble in standardized junction format
#' @examples
#' path <-  system.file("extdata", "", package = "splice2neo")
#' leafcutter_juncs <- leafcutter_transform(path)
#' leafcutter_juncs
#'
#' @import readr
#' @export
leafcutter_transform <- function(path) {

  # use strand information from bam.junc file as information in the cluster column in counts table seems to be wrong in some cases
  message("Importing Leafcutter files...")

  file.counts <- list.files(path, pattern = "perind.counts.gz")

  if(length(file.counts) == 0 ){stop("perind.counts.gz file is missing")}

  file.counts <- paste0(path, "/", file.counts)
  counts <- file.counts %>%
    import_leafcutter_counts() %>%
    transform_leafcutter_counts()%>%
    distinct(junc_id, .keep_all = TRUE)

  message("Importing Leafcutter files completed")

  return(counts)
}


