

# SPLADDER FUNCTION -------------------------------------------------------


#' Sorts columns of junction output file in the following order:
#' "junction_start", "junction_end", "strand", "chromosome", "Gene",
#' "class", "AS_event_ID", "junction_id"
#'
#' @param tib A tibble with the following as given in the description.
#'
#' @return A tibble with sorted columns as given above
#'
#' @examples
#'
#' unsorted_junc_df
#' sorted_junc_df <- sort_columns(unsorted_junc_df)
#' sorted_junc_df
#'
#'@import dplyr
#'
#'@export
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


#' Transforms events from alternative 3' splice sites from SPLADDER output
#' format into standardized junction format
#'
#' @param tib A tibble in SPLADDER output format
#'
#' @return A tibble in standardized junction format
#'
#' @examples
#' spladder_output.a3ss
#' transformed_spladder.a3ss <- spladder.transform.a3ss(spladder_output.a3ss)
#' transformed_spladder.a3ss
#'
#' @import dplyr
#' @export
spladder.transform.a3ss <- function(tib) {
  tib %>%
    mutate(
      junc_id1 = ifelse(
        .data$strand == "+",
        paste(.data$contig, .data$exon_const_end, .data$exon_alt1_start, .data$strand, sep = "_"),
        paste(.data$contig, .data$exon_alt1_end, .data$exon_const_start, .data$strand, sep = "_")
      ),
      junc_id2 = ifelse(
        .data$strand == "+",
        paste(.data$contig, .data$exon_const_end, .data$exon_alt2_start, .data$strand, sep = "_"),
        paste(.data$contig, .data$exon_alt2_end, .data$exon_const_start, .data$strand, sep = "_")
      ),
      junc_id = paste(.data$junc_id1, .data$junc_id2, sep = ";"),
      class = "A3SS",
      AS_event_ID = event_id
    ) %>%
    tidyr::separate_rows(.data$junc_id, sep = ";") %>%
    dplyr::select(.data$junc_id, .data$gene_name, .data$class, .data$AS_event_ID)
}


#' Transforms events from alternative 5' splice sites from SPLADDER output
#' format into standardized junction format
#'
#' @param tib A tibble in SPLADDER output format
#'
#' @return A tibble in standardized junction format
#'
#' @examples
#'spladder_output.a5ss
#'transformed_spladder.a5ss <- spladder.transform.a5ss(spladder_output.a5ss)
#'transformed_spladder.a5ss
#'
#'@import dplyr
#' @export
spladder.transform.a5ss <- function(tib) {
  tib %>%
    mutate(
      junc_id1 = ifelse(
        strand == "+",
        paste(contig, exon_alt1_end, exon_const_start, strand, sep = "_"),
        paste(contig, exon_const_end, exon_alt1_start, strand, sep = "_")
      ),
      junc_id2 = ifelse(
        strand == "+",
        paste(contig, exon_alt2_end, exon_const_start, strand, sep = "_"),
        paste(contig, exon_const_end, exon_alt2_start, strand, sep = "_")
      ),
      junc_id = paste(junc_id1, junc_id2, sep = ";"),
      class = "A5SS",
      AS_event_ID = event_id
    ) %>%
    tidyr::separate_rows(junc_id, sep = ";") %>%
    dplyr::select(junc_id, gene_name, class, AS_event_ID)
}

#' Transforms events resulting from exon skipping from SPLADDER output
#' format into standardized junction format
#'
#' @param tib A tibble in SPLADDER output format
#'
#' @return A tibble in standardized junction format
#'
#' @examples
#'spladder_output.exonskip
#'transformed_spladder.exonskip <- spladder.transform.exon_skipping(
#'  spladder_output.exonskip)
#'transformed_spladder.exonskip
#'
#'@import dplyr
#' @export
spladder.transform.exon_skipping <- function(tib) {
  tib %>%
    mutate(
      junc_id1 = paste(contig, exon_pre_end, exon_start, strand, sep = "_"),
      junc_id2 = paste(contig, exon_end, exon_aft_start, strand, sep = "_"),
      junc_id3 = paste(contig, exon_pre_end, exon_aft_start, strand, sep = "_"),
      junc_id = paste(junc_id1, junc_id2, junc_id3, sep = ";"),
      class = "cassette_exon",
      AS_event_ID = event_id
    ) %>%
    separate_rows(junc_id, sep = ";") %>%
    dplyr::select(junc_id, gene_name, class, AS_event_ID)
}


#' Transforms events resulting from intron retention from SPLADDER output
#' format into standardized junction format
#'
#' @param tib A tibble in SPLADDER output format
#'
#' @return A tibble in standardized junction format
#'
#' @examples
#'spladder_output.intronreten
#'transformed_spladder.intronreten <- spladder.transform.intron_retention(
#'  spladder_output.intronreten)
#'transformed_spladder.intronreten
#'
#'@import dplyr
#' @export
spladder.transform.intron_retention <- function(tib) {
  tib %>%
    mutate(
      junc_id1 = paste(contig, exon1_end, intron_start, strand, sep = "_"),
      junc_id2 = paste(contig, intron_end, exon2_start, strand, sep = "_"),
      class = "intron_retention",
      AS_event_ID = event_id,
      junc_id = paste(junc_id1, junc_id2, sep = ";")
    ) %>%
    separate_rows(junc_id, sep = ";") %>%
    dplyr::select(junc_id, gene_name, class, AS_event_ID)
}


#' Transforms events resulting from mutually exclusive exons from SPLADDER output
#' format into standardized junction format
#'
#' @param tib A tibble in SPLADDER output format
#'
#' @return A tibble in standardized junction format
#'
#' @examples
#'spladder_output.mutexon
#'transformed_spladder.mutexon <- spladder.transform.mutex_exon(
#'  spladder_output.mutexon)
#'transformed_spladder.mutexon
#'
#'@import dplyr
#' @export
spladder.transform.mutex_exon <- function(tib) {
  tib %>%
    mutate(
      junc_id1 = paste(contig, exon_pre_end, exon1_start, strand, sep = "_"),
      junc_id2 = paste(contig, exon1_end, exon_aft_start, strand, sep = "_"),
      junc_id3 = paste(contig, exon_pre_end, exon2_start, strand, sep = "_"),
      junc_id3 = paste(contig, exon2_end, exon_aft_start, strand, sep = "_"),
      junc_id = paste(junc_id1, junc_id2, junc_id3, sep = ";"),
      class = "mutex_exon",
      AS_event_ID = event_id
    ) %>%
    separate_rows(junc_id, sep = ";") %>%
    dplyr::select(junc_id, gene_name, class, AS_event_ID)
}


#' Transforms SPLADDER output into standardized junction format
#'
#' @param l A list with tibbles. Each tibble is a Spladder output for "A5SS",
#'  "A3SS", "cassette_exon", "intron_retention", "mutex_exons".
#'
#' @return A tibble in standardized junction format, combining all alternative
#'   splicing classes that are covered by Spladder
#'
#' @examples
#'spladder_output
#'transformed_spladder <- spladder.transform.format(
#'  spladder_output)
#'transformed_spladder
#'
#'@import dplyr
#' @export
spladder.transform.format <- function(l) {

  if(!all(c("A5SS","A3SS", "cassette_exon", "intron_retention",
            "mutex_exons") %in% names(l))){
    stop("The input list must contain the tibbles: A5SS, A3SS, cassette_exon,
    intron_retention, mutex_exons")
  }

  df <- bind_rows(
    spladder.transform.a5ss(l$A5SS),
    spladder.transform.a3ss(l$A3SS),
    spladder.transform.exon_skipping(l$cassette_exon),
    spladder.transform.intron_retention(l$intron_retention),
    spladder.transform.mutex_exon(l$mutex_exons)
  )
  df %>%
    tidyr::separate(
      junc_id,
      into = c("chromosome", "junction_start", "junction_end", "strand"),
      sep = "_",
      remove = F
    ) %>%
    mutate(
      junction_start = as.numeric(junction_start),
      junction_end = as.numeric(junction_end)
    ) %>%
    dplyr::rename(., Gene = gene_name) %>%
    sort_columns()

}


#' Imports SPLADDER output from a given path
#'
#' @param path The path to a folder with spladder output
#'
#' @return A list with tibbles. Each tibble is a Spladder output for "A5SS",
#'  "A3SS", "cassette_exon", "intron_retention", "mutex_exons".
#'
#'
#' @import readr
#' @export
import_spladder <- function(path){
  files <- list.files(path, "confirmed.txt.gz")
  path_files <- paste(path, files ,sep = "/" )
  files <- lapply(path_files, function(x) if(file.exists(x)) {read_delim(x, delim = "\t")})
  names(files) <- c("A3SS", "A5SS", "cassette_exon", "intron_retention", "mutex_exons")
  return(files)
}


#' Imports SPLADDER output from a given path and transforms it into standardized
#' junction format
#'
#' @param path The path to a folder with spladder output
#'
#' @return A tibble in standardized junction format, combining all alternative
#'   splicing classes that are covered by Spladder
#'
#'
#' @import readr
#' @export
spladder_transform <- function(path){
  dat <- import_spladder(path)
  juncs <- spladder.transform.format(dat)
  return(juncs)
}

# LEAFCUTTER FUNCTIONS -----------------------------------------------------

#' Imports "_perind.counts.gz" from Leafcutter outpu.
#'
#' @param file.name The path to the file
#'
#' @return A tibble with columns "intron_cluster" and "counts
#'
#'
#' @import readr
#' @export
import_leafcutter.counts <- function(file.name) {
  if(!file.exists(file.name)){
    stop("Aligned.out.bam.junc file is missing")
  }
  myData <- read_delim(file.name,
                       delim = " ",
                       skip = 1,
                       col_names = F)
  col.nams <- c("intron_cluster", "counts")
  colnames(myData) <- col.nams
  return(myData)
}


#' Transforms Leafcutter counts file
#'
#' @param tib Leafcutter counts file as tibble
#'
#' @return A tibble with columns "chromosome", "junction_start", "junction_end",
#'  "cluster"
#'
#'
#' @import dplyr
#' @export
transform_leafcutter.counts <- function(tib) {
  tib %>%
    separate(
      intron_cluster,
      sep = ":",
      into = c("chromosome", "junction_start", "junction_end", "cluster")
    ) %>%
    mutate(
      junction_start = as.numeric(junction_start),
      junction_end = as.numeric(junction_end),
      junc_id = paste(chromosome, junction_start, junction_end, sep = "_")
    )
}


#' Imports "_Aligned.out.bam.junc" file from Leafcutter output
#'
#' @param file.name The path to the file
#'
#' @return A tibble with columns "chromosome", "junction_start", "junction_end",
#' "dot", "count", "strand"
#'
#'
#' @import readr
#' @export
import_leafcutter.bam <- function(file.name) {
  if(!file.exists(file.name)){
    stop("Aligned.out.bam.junc file is missing")
  }
  bam.file <- read_delim(file.name,
                         delim = "\t", col_names = F)
  bam.cols <- c("chromosome",
                "junction_start",
                "junction_end",
                "dot",
                "count",
                "strand")
  colnames(bam.file) <- bam.cols
  return(bam.file)
}

#' Transforms Leafcutter bam junc file
#'
#' @param tib Leafcutter bam junc file as tibble
#'
#' @return A tibble including a column with the junction id
#'
#'
#' @import dplyr
#' @export
transform_leafcutter.bam <- function(tib) {
  tib %>%
    mutate(junc_id = paste(chromosome, junction_start, junction_end + 1,
                           sep = "_")) %>%
    select(junc_id, strand)
}


#' Transforms Leafcutter intermediate files into standardized format
#'
#' @param tib a tibble from leafcutter intermediate format
#'
#' @return A tibble in standardized junction format
#'
#'
#' @import dplyr
#' @export
leafcutter.generate.output <- function(tib) {
  tib <- tib %>%
    mutate(
      Gene = NA,
      class = NA,
      AS_event_ID = NA,
      junc_id = paste(chromosome, junction_start, junction_end, strand,
                          sep = "_")
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

#' Imports "_perind.counts.gz" from Leafcutter output and tranforms it
#' into standardized output format
#'
#' @param path The path to leafcutter output
#'
#' @return A tibble in standardized junction format
#'
#'
#' @import readr
#' @export
leafcutter_transform <- function(path) {
  file.counts <- list.files(path, pattern = "_perind.counts.gz")
  file.counts <- paste0(path, "/", file.counts)
  counts <- file.counts %>%
    import_leafcutter.counts() %>%
    transform_leafcutter.counts()
  file.bam <- list.files(path, pattern = "_Aligned.out.bam.junc")
  file.bam <- paste0(path, "/", file.bam)
  juncs <- file.bam %>%
    import_leafcutter.bam() %>%
    transform_leafcutter.bam()
  dat <- left_join(counts, juncs, by = "junc_id")
  dat_out <- leafcutter.generate.output(dat)
  return(dat_out)
}



# COMBINED DATASET --------------------------------------------------------


#' Combines tibbles with junctions from Spladder and Leafcutter
#'
#' @param leafcutter_juncs A tibble the junctions identified by Leafcutter in
#'   standardized format
#' @param spladder_juncs A tibble the junctions identified by Spladder in
#'   standardized format
#'
#' @return A combined table with unique junctions. The columns RNA_tool
#'   contains information which tools identified the given junction
#'
#' @examples
#' dat.combined <- generate_combined_dataset(spladder_juncs,
#'   leafcutter_juncs)
#' colnames(dat.combined)
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

#' This is wrapper function to map the information if a junction predicted from WES data was found in RNA-seq by leafcutter or spladder
#'
#' @param path_to_spladder The path to the results from RNA-seq analysis with spladder
#' @param path_to_leafcutter The path to the results from RNA-seq analysis with leafcutter
#' @param mutation_juncs The junction-transcript centric data predicted from WES data.
#'
#' @return The junction-transcript centric data predicted from WES data is extedended by the information if a respective aberrant junctions was
#' identified by spladder or leafcutter. (`identified_by_leafcutter`, `identified_by_spladder`)
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
