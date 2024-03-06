
library(splice2neo)
library(tidyverse)
# load genome of choice
library(BSgenome.Hsapiens.UCSC.hg19)
library(AnnotationDbi)

# # this is an example of a transcript database (see below for details)
# gtf_url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz"
# # parse GTF file as txdb object
# txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_url)
#
# transcripts <-
#   GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)
# transcripts_gr <- GenomicFeatures::transcripts(txdb, columns = c("gene_id", "tx_id", "tx_name"))
# # Build a GRangesList with cds composed of individual exon ranges
# cds <- GenomicFeatures::cdsBy(txdb, by = c("tx"), use.name = TRUE)

txdb_file <- "/projects/SUMMIT/WP1.1/splice/results/splice_analysis_v09.txdb.sqlite"
txdb <- loadDb(txdb_file)
transcripts <- GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)
exons_gr <- unlist(transcripts)
transcripts_gr <- GenomicFeatures::transcripts(txdb, columns = c("gene_id", "tx_id", "tx_name"))

# Build a GRangesList with cds composed of individual exon ranges
cds <- GenomicFeatures::cdsBy(txdb, by = c("tx"), use.name = TRUE)

gene_table <- read_delim("/projects/SUMMIT/WP1.1/alternative_splicing/data/canonical_junctions/genes_strand.tsv")


#gene_table <- read_delim(fs::path("test_data/", "gene_table.tsv"))

setwd("/projects/SUMMIT/WP1.1/alternative_splicing/splice2neo_new/")

canonical_juncs <- bed_to_junc(bed_file = "test_data/test_canonical.bed", type = "exon-exon")

dat_leafcutter <-
  leafcutter_transform(path = "test_data/leafcutter/")
# supported events: exon_skip, intron_retention, alt_3prime, alt_5prime, mutex_exons
dat_spladder <-
  spladder_transform(path = "test_data/spladder/")
dat_rna <-
  generate_combined_dataset(spladder_juncs = dat_spladder, leafcutter_juncs = dat_spladder)


dat_spliceai <-
  parse_spliceai(vcf_file = "test_data/spliceAI.vcf")
dat_splicai_formatted <- format_spliceai(dat_spliceai, gene_table = gene_table)
dat_spliceai_annotated <-
  annotate_mut_effect(effect_df = dat_splicai_formatted,
                      transcripts = transcripts,
                      transcripts_gr = transcripts_gr,
                      gene_mapping = TRUE)



dat_mmsplice <- parse_mmsplice(infile = "test_data/MMsplice.csv")
dat_mmsplice_annotated  <-
  annotate_mmsplice(mmsplice_df = dat_mmsplice, transcripts = transcripts)


dat_spliceai_annotated_unique <- unique_mut_junc(dat_spliceai_annotated)
dat_mmsplice_annotated_unique <- unique_junc_mmsplice(dat_mmsplice_annotated)


dat_mut <- combine_mut_junc(list(
  "spliceai" = dat_spliceai_annotated_unique,
  "mmsplice" = dat_mmsplice_annotated_unique
))

dat_mut <- dat_mut %>%
  mutate(is_canonical = is_canonical(junc_id, ref_junc = canonical_juncs, exons_gr = transcripts)) %>%
  mutate(is_in_rnaseq = is_in_rnaseq(junc_id, rna_juncs = dat_rna$junc_id))

# junction df --> contains not selected junction ids because we predict all potential junctions per mutation
dat_for_requantification <- dat_mut %>%
  filter(!is_canonical)

dat_for_requantification_cts <- dat_for_requantification %>%
  add_context_seq(size = 400, bsg = BSgenome.Hsapiens.UCSC.hg19, transcripts = transcripts)

dat_easyquant <- dat_for_requantification_cts %>%
  transform_for_requant()

dat_for_requantification_cts_peptide <-
  dat_for_requantification_cts  %>%
  add_peptide(flanking_size = 13, bsg = BSgenome.Hsapiens.UCSC.hg19, cds = cds)

dat_for_requantification_cts_peptide_req <- map_requant("test_data/easyquant/", dat_for_requantification_cts_peptide)


# TODO: wee need all potential transcripts for easyquant and not pre-selected ones

table(dat_for_requantification_cts_peptide$cts_id %in% easyq$name )

dat_for_requantification_cts_peptide %>% filter(! cts_id %in% easyq$name) %>%
  dplyr::select(junc_id, mut_id, tx_id, cts_id, is_canonical)


x <- dat_for_requantification_cts  %>% filter(!cts_id %in% easyq$name)
x$mut_id %in% dat2keep$mut_id

x %>% filter(mut_id == "chr10_91528535_G_A") %>% dplyr::select(junc_id, tx_id, spliceai_event_type, cts_id)
dat2keep %>% filter(mut_id == "chr10_91528535_G_A") %>% dplyr::select(junc_id, tx_id, spliceai_event_type)
