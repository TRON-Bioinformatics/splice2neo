
#' Sorts columns of junction output file in the following order: "junction_start", "junction_end", "strand", "chromosome", "Gene", "class", "AS_event_ID", "junction_id"
#'
#' @param tib A tibble with the following as given in the description.
#'
#' @return A tibble with sorted columns as given above
#'
#' @examples
#'
#' library(tibble)
#' tib <- tibble("strand" = "+", "chromosome" = "chr1", "Gene" = NA,"junction_start" = 0, "junction_end" = 111,  "class" = "intron_retention", "AS_event_ID" = NA, "junction_id" = "chr1_0_111_+")
#' tib <- sort_columns(x)
#' print(tib)
#'
#' @export
sort_columns <- function(tib){
  column_order <- c("junction_start", "junction_end", "strand", "chromosome",
                    "Gene", "class",
                    #"junction_info",
                    "AS_event_ID", "junction_id")
  tib[column_order]
}


#' Transforms events from alternative 3' splice sites from SPLADDER output format into standardized junction format
#'
#' @param tib A tibble in SPLADDER output format
#'
#' @return A tibble in standardized junction format
#'
#' @examples
#'
#' library(tibble)
#' x <- tibble("strand" = "+", "chromosome" = "chr1", "Gene" = NA,"junction_start" = 0, "junction_end" = 111,  "class" = "intron_retention", "AS_event_ID" = NA, "junction_id" = "chr1_0_111_+")
#' x1 <- sort_columns(x)
#' print(x1)
#'
#' @export
spladder.transform.a3ss <- function(tib){
  tib %>%
    mutate(junc_id1 = ifelse(strand == "+", paste(contig, exon_const_end,exon_alt1_start, strand, sep = "_"),
                             paste(contig, exon_alt1_end, exon_const_start, strand, sep = "_")),
           junc_id2 = ifelse(strand == "+", paste(contig, exon_const_end,exon_alt2_start, strand, sep = "_"),
                             paste(contig, exon_alt2_end, exon_const_start, strand, sep = "_")),
           junction_id = paste(junc_id1, junc_id2, sep = ";"),
           class = "A3SS",
           AS_event_ID = event_id) %>%
    separate_rows(junction_id, sep = ";") %>%
    dplyr::select(junction_id, gene_name, class, AS_event_ID)
}

#' Transforms events from alternative 5' splice sites from SPLADDER output format into standardized junction format
#'
#' @param x A tibble in SPLADDER output format
#'
#' @return A tibble in standardized junction format
#'
#' @examples
#'
#' library(tibble)
#' x <- tibble("strand" = "+", "chromosome" = "chr1", "Gene" = NA,"junction_start" = 0, "junction_end" = 111,  "class" = "intron_retention", "AS_event_ID" = NA, "junction_id" = "chr1_0_111_+")
#' x1 <- sort_columns(x)
#' print(x1)
#'
#' @export
spladder.transform.a5ss <- function(df){
  df %>%
    mutate(junc_id1 = ifelse(strand == "+", paste(contig, exon_alt1_end, exon_const_start, strand, sep = "_"),
                             paste(contig, exon_const_end, exon_alt1_start, strand, sep = "_")),
           junc_id2 = ifelse(strand == "+", paste(contig, exon_alt2_end, exon_const_start, strand, sep = "_"),
                             paste(contig, exon_const_end, exon_alt2_start, strand, sep = "_")),
           junction_id = paste(junc_id1, junc_id2, sep = ";"),
           class = "A5SS",
           AS_event_ID = event_id) %>%
    separate_rows(junction_id, sep = ";") %>%
    dplyr::select(junction_id, gene_name, class, AS_event_ID)
}

# cassette_exon / exon skipping
spladder.transform.exon_skipping <- function(df){
  df %>%
    mutate(junc_id1 = paste(contig, exon_pre_end, exon_start, strand, sep = "_"),
           junc_id2 = paste(contig, exon_end, exon_aft_start, strand, sep = "_"),
           junc_id3 = paste(contig, exon_pre_end, exon_aft_start, strand, sep = "_"),
           junction_id = paste(junc_id1, junc_id2, junc_id3, sep = ";"),
           class = "cassette_exon",
           AS_event_ID = event_id) %>%
    separate_rows(junction_id, sep = ";") %>%
    dplyr::select(junction_id, gene_name, class, AS_event_ID)
}


# intron retention
# downstream exon-intron boundary
spladder.transform.intron_retention <- function(df){
  df %>%
    mutate(junc_id1 = paste(contig, exon1_end, intron_start, strand, sep = "_"),
           junc_id2 = paste(contig, intron_end, exon2_start, strand, sep = "_"),
           class = "intron_retention",
           AS_event_ID = event_id,
           junction_id = paste(junc_id1, junc_id2, sep = ";")) %>%
    separate_rows(junction_id, sep = ";") %>%
    dplyr::select(junction_id, gene_name, class, AS_event_ID)
}

# mutually exclusive exons
spladder.transform.mutex_exon <- function(df){
  df %>%
    mutate(junc_id1 = paste(contig, exon_pre_end, exon1_start, strand, sep = "_"),
           junc_id2 = paste(contig, exon1_end, exon_aft_start, strand, sep = "_"),
           junc_id3 = paste(contig, exon_pre_end, exon2_start, strand, sep = "_"),
           junc_id3 = paste(contig, exon2_end, exon_aft_start, strand, sep = "_"),
           junction_id = paste(junc_id1, junc_id2, junc_id3, sep = ";"),
           class = "mutex_exon",
           AS_event_ID = event_id) %>%
    separate_rows(junction_id, sep = ";") %>%
    dplyr::select(junction_id, gene_name, class, AS_event_ID)
}


spladder.transform.format <- function(list.dfs){
  df <- bind_rows(
    spladder.transform.a5ss(list.dfs$A5SS),
    spladder.transform.a3ss(list.dfs$A3SS),
    spladder.transform.exon_skipping(list.dfs$cassette_exon),
    spladder.transform.intron_retention(list.dfs$intron_retention),
    spladder.transform.mutex_exon(list.dfs$mutex_exons)
  )
  df %>%
    separate(junction_id, into = c("chromosome", "junction_start", "junction_end", "strand"), sep = "_", remove = F) %>%
    mutate(junction_start = as.numeric(junction_start), junction_end = as.numeric(junction_end))%>%
    dplyr::rename(., Gene = gene_name ) %>%
    sort_columns()

}


# LEAFCUTTER --------------------------------------------------------------

# import filtered leafcutter output
leafcutter.import <- function(file.name){
  myData <- read_delim(file.name, delim = " ", skip = 1,
                       col_names = F)
  col.nams <- c("intron_cluster", "counts")
  colnames(myData) <- col.nams
  return(myData)
}

# import leafcutter junction file
import.bam <- function(file.name){
  bam.file <- read_delim(file.name,
                         delim = "\t", col_names = F)
  bam.cols <- c("chromosome", "junction_start", "junction_end",
                "dot", "count", "strand")
  colnames(bam.file) <- bam.cols
  return(bam.file)
}

# adds junc id for matching with strand information to filtered leafcutter output
leafcutter.transform.input <- function(x){
  x %>%
    separate(intron_cluster, sep = ":", into = c("chromosome","junction_start","junction_end", "cluster")) %>%
    mutate(junction_start = as.numeric(junction_start),
           junction_end = as.numeric(junction_end),
           junc_id = paste(chromosome, junction_start, junction_end, sep = "_"))
}

# adds junc id for matching with strand information to leafcutter junctions output
transform.bam <- function(x){
  x %>%
    mutate(junc_id = paste(chromosome, junction_start, junction_end + 1,
                           sep = "_")) %>%
    select(junc_id, strand)
}

# transforms into standardised output
leafcutter.generate.output <- function(x){
  x <- x %>%
    mutate(Gene = NA, class =NA, AS_event_ID = NA,
           junction_id = paste(chromosome, junction_start, junction_end, strand, sep ="_")) %>%
    select(junction_start,junction_end,strand,chromosome,Gene,class,AS_event_ID,junction_id)
  return(x)
}

# combines junction/events tibbles of replicates into one tibble and groups by patient
combine_reps2tibble <- function(list.tibbles.replicates){
  list.tibbles.replicates %>%
    bind_rows(.id = "replicate") %>%
    mutate(tumour = gsub("_R[1,2]", "", replicate), junc_tumor = paste0(junction_id, "&",tumour)) %>%
    group_by(junc_tumor) %>%
    summarise(replicates = paste(replicate, collapse=","), junction_start = junction_start, junction_end = junction_end,
              strand = strand, Gene = Gene, class = class, AS_event_ID = AS_event_ID, junction_id = junction_id,
              tumour = tumour) %>%
    separate(replicates, into = c("replicate1", "replicate2"), sep = ",") %>%
    mutate(replicate1 = ifelse(is.na(replicate1),F,T),replicate2 = ifelse(is.na(replicate2),F,T)) %>%
    ungroup() %>%
    select(-junc_tumor) %>%
    distinct() %>%
    group_by(tumour)
}
