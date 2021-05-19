## code to prepare `toy_transcript` dataset goes here
## This dataset contains human genomic annotations for the region hg19:
##
##  chr2  152000000 - 180000000 # example data for SpliceAI
##  chr17 41100000 - 41280000   # example data for MMSplice

# 17:41197805
# 17:41277297

require(tidyverse)

gtf_url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz"

txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_url)

# filter for a subset of chromosomes. See https://bioconductor.org/packages/devel/bioc/vignettes/GenomicFeatures/inst/doc/GenomicFeatures.pdf
# head(GenomeInfoDb::seqlevels(txdb))

GenomeInfoDb::seqlevels(txdb) <- c("chr2")

# GenomicFeatures::seqlevels(txdb) <- GenomicFeatures::seqlevels0(txdb)

# Build a GRangesList with transcripts composed of exon ranges
transcripts <- GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)
# pryr::object_size(transcripts)

# Build GRanges of whole transcripts
transcripts_gr <- GenomicFeatures::transcripts(txdb)

# Build GRanges of all exons
exons_gr <- GenomicFeatures::exons(txdb, columns = "exon_name", use.names = TRUE)


# Restict toy data to hg19 chr2 152000000 - 180000000
sub_gr <- GenomicRanges::GRanges(c(
  "chr2:152000000-180000000",
  "chr17:41100000-41280000"
  ))

toy_transcripts <- transcripts %>%
  IRanges::subsetByOverlaps(sub_gr)

toy_transcripts_gr <- transcripts_gr %>%
  IRanges::subsetByOverlaps(sub_gr)

toy_exons_gr <- exons_gr %>%
  IRanges::subsetByOverlaps(sub_gr)

# Add data to package
# usethis::use_data(txdb, overwrite = TRUE)
usethis::use_data(toy_transcripts, overwrite = TRUE)
usethis::use_data(toy_transcripts_gr, overwrite = TRUE)
# usethis::use_data(toy_exons_gr, overwrite = TRUE)



