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
    mutate(
      junction_start = as.numeric(junction_start),
      junction_end = as.numeric(junction_end),
      fake_junc_id = paste(chromosome, junction_start, junction_end, sep = "_")
    )
}


#' Imports "_Aligned.out.bam.junc" file from LeafCutter output
#'
#' @param file.name The path to the file
#'
#' @return A tibble with columns "chromosome", "junction_start", "junction_end",
#' "dot", "count", "strand"
#'
#'
#' @import readr
#' @export
import_leafcutter_bam <- function(file.name) {
  if(!file.exists(file.name)){
    stop("Aligned.out.bam.junc file is missing")
  }
  bam.file <- read_delim(file.name,
                         delim = "\t",
                         col_names = F,
                         col_types = cols(.default = "c"))
  bam.cols <- c("chrom",
                "chromStart",
                "chromEnd",
                "name",
                "score",
                "strand",
                "thickStart",
                "thickEnd",
                "itemRgb",
                "blockCount",
                "blockSizes",
                "blockStarts")
  colnames(bam.file) <- bam.cols
  return(bam.file)
}

#' Transforms LeafCutter bam junc file
#'
#' @param tib LeafCutter bam junc file as tibble
#'
#' @return A tibble including a column with the junction id
#'
#'
#' @import dplyr
#' @export
transform_leafcutter_bam <- function(tib) {
  tib %>%
    separate(blockSizes, into = c("start_off", "end_off"), sep = ",",
             remove = FALSE) %>%
    mutate(exact_start = as.numeric(chromStart) + as.numeric(start_off),
           exact_end = as.numeric(chromEnd) - as.numeric(end_off) + 1 )%>%
    mutate(fake_junc_id = paste(chrom, exact_start, exact_end , sep = "_")) %>%
    select(fake_junc_id, strand)
}


#' Transforms LeafCutter intermediate files into standardized format
#'
#' @param tib a tibble from LeafCutter intermediate format
#'
#' @return A tibble in standardized junction format
#'
#'
#' @import dplyr
#' @export
leafcutter_transform_format <- function(tib) {
  tib <- tib %>%
    mutate(
      Gene = NA,
      class = NA,
      AS_event_ID = NA,
      junc_id = generate_junction_id(chromosome, junction_start, junction_end, strand)
    ) %>%
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
  return(tib)
}

#' Imports "*perind.counts.gz" and "*bam.junc" from LeafCutter output and transforms the raw output
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
    transform_leafcutter_counts()

  file.bam <- list.files(path, pattern = "bam.junc")
  if(length(file.bam) == 0 ){stop("bam.junc file is missing")}
  if(length(file.bam) > 2 ){stop("There is more than one bam.junc file present in the folder")}

  file.bam <- paste0(path, "/", file.bam)
  juncs <- file.bam %>%
    import_leafcutter_bam() %>%
    transform_leafcutter_bam()
  dat <- left_join(counts, juncs, by = "fake_junc_id")
  dat_out <- leafcutter_transform_format(dat)

  message("Importing Leafcutter files completed")

  return(dat_out)
}


