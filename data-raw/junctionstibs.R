## code to prepare junction example datasets goes here

library(tibble)


unsorted_junc_df = tibble(
  "strand" = "+", "chromosome" = "chr1", "Gene" = NA,
  "junction_start" = 0, "junction_end" = 111,  "class" = "intron_retention",
  "AS_event_ID" = NA, "junc_id" = "chr1_0_111_+")





#
# spladder_output.a5ss <- tibble(
#   "chrm" = c("chr1", "chr1"),
#   "strand" = c("+", "-"),
#   "event_id" = c("alt_5prime.1","alt_5prime.1825"),
#   "is_annotated" = c("3","3"),
#   "gene_name" = c("ENSG00000171163.16_6",""),
#   "e1_start" = c("763186","934906"),
#   "e1_end" = c("763229","934993"),
#   "e2_start" = c("763230","935072"),
#   "e2_end" = c("763233","935245"),
#   "e3_start" = c("764383","935246"),
#   "e3_end" = c("764484","935361"),
#   "sample_Aligned.out_sorted:exon_diff_cov" = 0,
#   "sample_Aligned.out_sorted:exon_const_cov" = 0,
#   "sample_Aligned.out_sorted:intron1_conf" = 0,
#   "sample_Aligned.out_sorted:intron2_conf" = 0,
#   "sample_Aligned.out_sorted:psi" = NA
# )
#
#
#
# spladder_output.exonskip <- tibble(
#   "chrm" = "chr1",
#   "strand" = "+",
#   "event_id" = "alt_5prime.1",
#   "gene_name" = "ENSG00000171163.16_6",
#   "e1_start" = "59228292",
#   "e1_end" = "59228349",
#   "e2_start" = "59228549",
#   "e2_end" = "59228749",
#   "e3_start" = "59229049",
#   "e3_end" = "59229949",
#   "sample_Aligned.out_sorted:exon_diff_cov" = 0,
#   "sample_Aligned.out_sorted:exon_const_cov" = 0,
#   "sample_Aligned.out_sorted:intron1_conf" = 0,
#   "sample_Aligned.out_sorted:intron2_conf" = 0,
#   "sample_Aligned.out_sorted:psi" = NA
# )
#
#
# spladder_output.intronreten <- tibble(
#   "chrm" = "chr1",
#   "strand" = "-",
#   "event_id" = "alt_5prime.1",
#   "gene_name" = "ENSG00000171163.16_6",
#   "e1_start" = "244974",
#   "e1_end" = "249445",
#   "e2_start" = "249446",
#   "e2_end" = "249512",
#   "e3_start" = "249513",
#   "e3_end" = "249631",
#   "sample_Aligned.out_sorted:exon_diff_cov" = 4,
#   "sample_Aligned.out_sorted:exon_const_cov" = 1,
#   "sample_Aligned.out_sorted:intron1_conf" = 0,
#   "sample_Aligned.out_sorted:intron2_conf" = 0,
#   "sample_Aligned.out_sorted:psi" = NA
# )
#
#
# spladder_output.mutexon <- tibble(
#   "chrm" = "chr1",
#   "strand" = "-",
#   "event_id" = "alt_5prime.1",
#   "gene_name" = "ENSG00000171163.16_6",
#   "e1_start" = "244974",
#   "e1_end" = "249445",
#   "e2_start" = "249446",
#   "e2_end" = "249512",
#   "e3_start" = "249513",
#   "e3_end" = "249631",
#   "e4_start" = "249513",
#   "e4_end" = "249631",
#   "sample_Aligned.out_sorted:exon_diff_cov" = 4,
#   "sample_Aligned.out_sorted:exon_const_cov" = 1,
#   "sample_Aligned.out_sorted:intron1_conf" = 0,
#   "sample_Aligned.out_sorted:intron2_conf" = 0,
#   "sample_Aligned.out_sorted:psi" = NA
# )
#
#
# spladder_output <- list(
#   "A3SS" = spladder_output.a3ss,
#   "A5SS" = spladder_output.a5ss,
#   "cassette_exon" = spladder_output.exonskip,
#   "intron_retention" = spladder_output.intronreten,
#   "mutex_exons" = spladder_output.mutexon
# )
#
#
# leafcutter_bam_junc <- tibble(
#   chrom = c("chr1"),
#   chromStart= c("21877792", "59228253"),
#   chromEnd= c("21878258", "59228633"),
#   name= c("JUNC00000001", "JUNC00000002"),
#   score= c("8", "15"),
#   strand= c("-", "+"),
#   thickStart= c("21877792", "59228253"),
#   thickEnd= c("21878258", "59228633"),
#   itemRg= c("255,0,0", "255,0,0"),
#   blockCount= c("2", "2"),
#   blockSizes= c("98,94", "96,85"),
#   blockStarts= c("0,238", "0,196")
# )
#
# leafcutter_counts <- tibble(
#   intron_cluster  = c("chr1:21877890:21878165:clu_1_-", "chr1:59228349:59228549:clu_2_+"),
#   counts= c("8/421", "413/421")
# )


leafcutter_counts  <- readr::read_delim("inst/extdata/test_perind.counts.gz",  col_types = cols(.default = "c"))
leafcutter_bam_junc <- readr::read_delim("inst/extdata/test_Aligned.out.sorted.bam.junc", col_types = cols(.default = "c"),
                                         e)
spladder_output <- import_spladder("inst/extdata/")

# leafcutter_juncs <- feather::read_feather("inst/extdata/leafcutter_example.feather")
# spladder_juncs <- feather::read_feather("inst/extdata/leafcutter_example.feather")
canonical_juncs <- feather::read_feather("data-raw/canonical_example.feather")


usethis::use_data(unsorted_junc_df, overwrite = TRUE)
# usethis::use_data(spladder_output.a3ss, overwrite = TRUE)
# usethis::use_data(spladder_output.a5ss, overwrite = TRUE)
# usethis::use_data(spladder_output.exonskip, overwrite = TRUE)
# usethis::use_data(spladder_output.intronreten, overwrite = TRUE)
# usethis::use_data(spladder_output.mutexon, overwrite = TRUE)
usethis::use_data(spladder_output, overwrite = TRUE)
# usethis::use_data(leafcutter_juncs, overwrite = TRUE)
# usethis::use_data(spladder_juncs, overwrite = TRUE)
usethis::use_data(canonical_juncs, overwrite = TRUE)
usethis::use_data(leafcutter_counts, overwrite = TRUE)
usethis::use_data(leafcutter_bam_junc, overwrite = TRUE)
