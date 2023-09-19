# REGTOOLS FUNCTIONS -----------------------------------------------------

#' Imports Regtools junctions annotate table
#'
#' Ignores junctions without strand annotation.
#'
#' @param file.name The path to the file
#'
#' @return A tibble with regtools identified junctions
#'
#'
#' @import readr
#' @export
import_regtools_junc <- function(file.name) {
    if(!file.exists(file.name)){
            stop("Regtools junctions annotate output file is missing")
    }
    junc.file <- read_tsv(file.name,
                            show_col_types = FALSE,
                            col_types = cols(.default = "c"))
    # 20211222 (JoHa): This filtering might be required since in rare cases the strand information
    # might be ambiguous ("?" for unstranded (XS-tagged) RNA-seq data)
    junc.file <- junc.file %>%
        filter(strand %in% c('+', '-'))
    return(junc.file)
}

#' Transforms Regtools intermediate files into standardized format
#'
#' @param tib A tibble from Regtools intermediate format
#'
#' @return A tibble in standardized junction format
#'
#'
#' @import dplyr tidyr
#' @export
transform_regtools_junc <- function(tib) {
    # 20211222 (JoHa): Event classification based on regtools documentation (add link)
    names_events <- c('DA'  = 'canonical_junction',
                      'NDA' = 'cassette_exon',
                      'D'   = 'A3SS',
                      'A'   = 'A5SS',
                      'N'   = 'novel_acceptor_novel_donor')
    tib <- tib %>%
        mutate(class = names_events[anchor],
               AS_event_ID = NA,
               junc_id = generate_junction_id(chrom, start, end, strand),
               score = as.numeric(score)) %>%
        dplyr::select(start,
                      end,
                      strand,
                      chrom,
                      gene_ids,
                      class,
                      AS_event_ID,
                      junc_id,
                      score)

    junc.cols <- c("junction_start",
                   "junction_end",
                   "strand",
                   "chromosome",
                   "Gene",
                   "class",
                   "AS_event_ID",
                   "junc_id",
                   "number_supporting_reads")
    colnames(tib) <- junc.cols
    return(tib)
}

#' Imports Regtools junctions annotate output and transforms the raw output
#' into standardized junction output format
#'
#'  - GitHub: https://github.com/griffithlab/regtools
#'  - Paper: https://doi.org/10.1038/s41467-023-37266-6
#'
#' Ignores junctions without strand annotation.
#'
#' @param path The path to the Regtools output file
#'
#' @return A tibble in standardized junction format
#'
#'
#' @import readr
#' @export
regtools_transform <- function(path) {
  file.junc <- path
  juncs     <- file.junc %>%
    import_regtools_junc() %>%
    transform_regtools_junc()
  return(juncs)
}
