## code to prepare `toy_genome` dataset goes here
## This dataset contains human genomic annotations for the region hg19
##  chr2  152000000 - 180000000 # example data for SpliceAI
##  chr17 41100000 - 41280000   # example data for MMSplice


require(tidyverse)
require(biomaRt)

# TODO: create toy genome based on download from embl --> not possible because of firewall
# ensembl <- useMart("ensembl")
# datasets <- listDatasets(ensembl)
# head(datasets)
# GENES = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "useast.ensembl.org")
#
# YouRange<- getBM(attributes = c("chromosome_name", "start_position", "end_position"),
#                  filters = c("chromosome_name","start","end"),
#                  values = list(chromosome_name = "4", start = 50000-3000,
#                                end = 50000 + 3000), mart = GENES)

path.to.genome <- "/projects/SUMMIT/WP1.1/alternative_splicing/annotations/GRCh37.primary_assembly.genome.fa"
# genome sequence

genome.seqs <- Biostrings::readDNAStringSet(path.to.genome)
names(genome.seqs) <- str_split_fixed(names(genome.seqs), " ", n = 2)[,1]

# select chromosome
#chroms <- which(names(genome.seqs) %in% c("chr2", "chr17"))
chroms <- which(names(genome.seqs) == "chr17")
genome <- genome.seqs[chroms]
genome <- Biostrings::DNAStringSet(genome[[1]][0:41300000])
names(genome) <- "chr17"


usethis::use_data(genome, overwrite = TRUE)
