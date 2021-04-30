## code to prepare junction example datasets goes here

library(tibble)


unsorted_junc_df = tibble(
  "strand" = "+", "chromosome" = "chr1", "Gene" = NA,
  "junction_start" = 0, "junction_end" = 111,  "class" = "intron_retention",
  "AS_event_ID" = NA, "junction_id" = "chr1_0_111_+")



spladder_output.a3ss <- tibble(
  "contig" = "chr1",
  "strand" = "-",
  "event_id" = "alt_3prime.1",
  "gene_name" = "ENSG00000171163.16_6",
  "exon_const_start" = "21878165",
  "exon_const_end" = "21878359",
  "exon_alt1_start" = "21877709",
  "exon_alt1_end" = "21877890",
  "exon_alt2_start" = "21877709",
  "exon_alt2_end" = "21877997",
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:exon_diff_cov" = 0,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:exon_const_cov" = 0,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:intron1_conf" = 0,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:intron2_conf" = 0,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:psi" = NA
)


spladder_output.a5ss <- tibble(
  "contig" = "chr1",
  "strand" = "-",
  "event_id" = "alt_5prime.1",
  "gene_name" = "ENSG00000171163.16_6",
  "exon_const_start" = "59228292",
  "exon_const_end" = "59228349",
  "exon_alt1_start" = "59222127",
  "exon_alt1_end" = "59222216",
  "exon_alt2_start" = "59222127",
  "exon_alt2_end" = "59222281",
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:exon_diff_cov" = 0,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:exon_const_cov" = 0,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:intron1_conf" = 0,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:intron2_conf" = 0,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:psi" = NA
)



spladder_output.exonskip <- tibble(
  "contig" = "chr1",
  "strand" = "+",
  "event_id" = "alt_5prime.1",
  "gene_name" = "ENSG00000171163.16_6",
  "exon_pre_start" = "59228292",
  "exon_pre_end" = "59228349",
  "exon_start" = "59228549",
  "exon_end" = "59228749",
  "exon_aft_start" = "59229049",
  "exon_aft_end" = "59229949",
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:exon_diff_cov" = 0,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:exon_const_cov" = 0,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:intron1_conf" = 0,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:intron2_conf" = 0,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:psi" = NA
)


spladder_output.intronreten <- tibble(
  "contig" = "chr1",
  "strand" = "-",
  "event_id" = "alt_5prime.1",
  "gene_name" = "ENSG00000171163.16_6",
  "exon1_start" = "244974",
  "exon1_end" = "249445",
  "intron_start" = "249446",
  "intron_end" = "249512",
  "exon2_start" = "249513",
  "exon2_end" = "249631",
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:exon_diff_cov" = 4,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:exon_const_cov" = 1,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:intron1_conf" = 0,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:intron2_conf" = 0,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:psi" = NA
)


spladder_output.mutexon <- tibble(
  "contig" = "chr1",
  "strand" = "-",
  "event_id" = "alt_5prime.1",
  "gene_name" = "ENSG00000171163.16_6",
  "exon_pre_start" = "244974",
  "exon_pre_end" = "249445",
  "exon1_start" = "249446",
  "exon1_end" = "249512",
  "exon2_start" = "249513",
  "exon2_end" = "249631",
  "exon_aft_start" = "249513",
  "exon_aft_end" = "249631",
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:exon_diff_cov" = 4,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:exon_const_cov" = 1,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:intron1_conf" = 0,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:intron2_conf" = 0,
  "IE000073B2_FF_capture_R1_Aligned.out_sorted:psi" = NA
)


spladder_output <- list(
  "A3SS" = spladder_output.a3ss,
  "A5SS" = spladder_output.a5ss,
  "cassette_exon" = spladder_output.exonskip,
  "intron_retention" = spladder_output.intronreten,
  "mutex_exons" = spladder_output.mutexon
)

usethis::use_data(unsorted_junc_df, overwrite = TRUE)
usethis::use_data(spladder_output.a3ss, overwrite = TRUE)
usethis::use_data(spladder_output.a5ss, overwrite = TRUE)
usethis::use_data(spladder_output.exonskip, overwrite = TRUE)
usethis::use_data(spladder_output.intronreten, overwrite = TRUE)
usethis::use_data(spladder_output.mutexon, overwrite = TRUE)
usethis::use_data(spladder_output, overwrite = TRUE)
