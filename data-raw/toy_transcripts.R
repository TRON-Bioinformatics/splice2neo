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

GenomeInfoDb::seqlevels(txdb) <- c("chr2", "chr17")

# GenomicFeatures::seqlevels(txdb) <- GenomicFeatures::seqlevels0(txdb)

# Build a GRangesList with transcripts composed of exon ranges
transcripts <- GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)
# pryr::object_size(transcripts)

# Build GRanges of whole transcripts
transcripts_gr <- GenomicFeatures::transcripts(txdb)

# # Build GRanges of all exons
# exons_gr <- GenomicFeatures::exons(txdb, columns = "exon_name", use.names = TRUE)

cds <- GenomicFeatures::cdsBy(txdb, by = c("tx"), use.name = TRUE)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Restict toy data to hg19 chr2 152000000 - 180000000 and chr17:41100000-41280000
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sub_gr <- GenomicRanges::GRanges(c(
  "chr2:152000000-180000000",
  "chr17:41100000-41280000"
))

toy_transcripts <- transcripts %>%
  IRanges::subsetByOverlaps(sub_gr)

toy_transcripts_gr <- transcripts_gr %>%
  IRanges::subsetByOverlaps(sub_gr)

toy_cds <- cds %>%
  IRanges::subsetByOverlaps(sub_gr)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# convert transcript name ENST ids by removing the version number suffix -------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
names(toy_transcripts) <- names(toy_transcripts) %>%
  str_extract("^ENST\\d{11}")
toy_transcripts <- toy_transcripts %>%
  S4Vectors::endoapply(function(x){x$exon_name = x$exon_name %>% str_extract("^ENSE\\d{11}"); x})

toy_transcripts_gr$tx_name <- toy_transcripts_gr$tx_name %>%
  str_extract("^ENST\\d{11}")

names(toy_cds) <- names(toy_cds) %>%
  str_extract("^ENST\\d{11}")


# Add data to package
usethis::use_data(toy_transcripts, overwrite = TRUE)
usethis::use_data(toy_transcripts_gr, overwrite = TRUE)
usethis::use_data(toy_cds, overwrite = TRUE)


################################################################################
## Toy junction IDs
################################################################################

# The toy junction ids were generated from the spliceAI example file with
# splice2neo v 0.0.1 as follows
#
#
# spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")
# df_raw <- parse_spliceai(spliceai_file)
# df <- format_spliceai(df_raw)
# annot_df <- annotate_spliceai_junction(df, toy_transcripts, toy_transcripts_gr)
# annot_df %>% select(junc_id) %>% distinct() %>% pull() %>% dput()
#
# toy_junc_id <- c("chr2_152389996_152392205_-", "chr2_152389996_152390729_-",
#   "chr2_152389955_152389956_-", "chr2_152388410_152392205_-", "chr2_152388410_152390729_-",
#   "chr2_179415981_179416357_-", "chr2_179415987_179415988_-", "chr2_179415000_179416357_-",
#   "chr2_179445336_179446207_-", "chr2_179446225_179446226_-", "chr2_179445336_179446633_-",
#   "chr2_179642187_179644565_+", "chr2_179642044_179642187_-", "chr2_179638062_179644565_+",
#   "chr2_179642146_179642147_-", "chr2_179642044_179642431_-")
#
# # Toy ENST IDs are unique ENST ID per junc ID that is affected by the junction
# toy_junc_id_enst <- c("ENST00000409198", "ENST00000409198", "ENST00000409198", "ENST00000409198",
#                       "ENST00000409198", "ENST00000342992", "ENST00000342992", "ENST00000342992",
#                       "ENST00000342992", "ENST00000342992", "ENST00000342992", "ENST00000342992",
#                       "ENST00000342992", "ENST00000342992", "ENST00000342992", "ENST00000342992")

toy_junc_id <- c("chr2:152389996-152392205:-", "chr2:152389996-152390729:-",
                 "chr2:152389955-152389956:-", "chr2:152388410-152392205:-", "chr2:152388410-152390729:-",
                 "chr2:179415981-179416357:-", "chr2:179415987-179415988:-", "chr2:179415000-179416357:-",
                 "chr2:179445336-179446207:-", "chr2:179446225-179446226:-", "chr2:179445336-179446633:-",
                 "chr2:179642044-179642187:-", "chr2:179642146-179642147:-", "chr2:179642044-179642431:-",
                 "chr2:152226533-152226534:+", "chr2:152222731-152222732:+", "chr2:152388410-152388411:-"
)

# Toy ENST IDs are unique ENST ID per junc ID that is affected by the junction
toy_junc_id_enst <- c("ENST00000409198", "ENST00000409198", "ENST00000397345", "ENST00000409198",
                      "ENST00000409198", "ENST00000342992", "ENST00000342992", "ENST00000342992",
                      "ENST00000342992", "ENST00000342992", "ENST00000342992", "ENST00000342992",
                      "ENST00000342992", "ENST00000342992", "ENST00000460812", "ENST00000460812", "ENST00000397345")



toy_junc_df <- tibble(
  junc_id = toy_junc_id,
  tx_id = toy_junc_id_enst
)

usethis::use_data(toy_junc_id, overwrite = TRUE)
usethis::use_data(toy_junc_id_enst, overwrite = TRUE)
usethis::use_data(toy_junc_df, overwrite = TRUE)


