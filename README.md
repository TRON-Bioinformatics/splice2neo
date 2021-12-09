
<!-- README.md is generated from README.Rmd. Please edit that file -->

# splice2neo

<!-- badges: start -->

[![R-CMD-check](https://github.com/TRON-Bioinformatics/splice2neo/workflows/R-CMD-check/badge.svg)](https://github.com/TRON-Bioinformatics/splice2neo/actions)
[![Codecov test
coverage](https://codecov.io/gh/TRON-Bioinformatics/splice2neo/branch/master/graph/badge.svg)](https://codecov.io/gh/TRON-Bioinformatics/splice2neo?branch=master)
[![](https://img.shields.io/badge/devel%20version-0.1.3-blue.svg)](https://github.com/TRON-Bioinformatics/splice2neo)
[![](https://img.shields.io/badge/lifecycle-experimental-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![](https://img.shields.io/github/last-commit/TRON-Bioinformatics/splice2neo.svg)](https://github.com/TRON-Bioinformatics/splice2neo/commits/dev)
<!-- badges: end -->

This package provides functions for the analysis of alternative splicing
junctions and their association with somatic mutations. It integrates
the output of several tools which predict splicing effects from
mutations or which detect expressed splice junctions from RNA-seq data
into a standardized splice junction format based on genomic coordinates.
Detected splice junctions can be filtered against canonical ones and
annotated with affected transcript sequences, CDS, and resulting peptide
sequences. The resulting tumor-specific splice junctions can encode
neoantigens.

Website: <https://tron-bioinformatics.github.io/splice2neo/>

## Installation

This R package is not yet on [CRAN](https://CRAN.R-project.org) or
[Bioconductor](https://www.bioconductor.org/). Therefore, you need to
install it from this repository.

``` r
## install.packages("remotes")
remotes::install_github("TRON-Bioinformatics/splice2neo")
```

## Example

This is a basic example of how to use some functions.

``` r
library(splice2neo)

# load human genome reference sequence
requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
```

### Example data

We start with some example splice junctions provided with the package.

``` r
junc_df <- dplyr::tibble(
  junc_id = toy_junc_id[c(1, 6, 10)]
)

junc_df
#> # A tibble: 3 × 1
#>   junc_id                   
#>   <chr>                     
#> 1 chr2_152389996_152392205_-
#> 2 chr2_179415981_179416357_-
#> 3 chr2_179446225_179446226_-
```

### Add transcripts

Next, we find the transcripts which are affected by the splice
junctions.

``` r
junc_df <- junc_df %>% 
  add_tx(toy_transcripts)

junc_df
#> # A tibble: 21 × 3
#>    junc_id                    tx_id           tx_lst      
#>    <chr>                      <chr>           <named list>
#>  1 chr2_152389996_152392205_- ENST00000409198 <GRanges>   
#>  2 chr2_152389996_152392205_- ENST00000172853 <GRanges>   
#>  3 chr2_152389996_152392205_- ENST00000397345 <GRanges>   
#>  4 chr2_152389996_152392205_- ENST00000427231 <GRanges>   
#>  5 chr2_152389996_152392205_- ENST00000618972 <GRanges>   
#>  6 chr2_152389996_152392205_- ENST00000413693 <GRanges>   
#>  7 chr2_152389996_152392205_- ENST00000603639 <GRanges>   
#>  8 chr2_152389996_152392205_- ENST00000604864 <GRanges>   
#>  9 chr2_152389996_152392205_- ENST00000420924 <GRanges>   
#> 10 chr2_179415981_179416357_- ENST00000342992 <GRanges>   
#> # … with 11 more rows
```

### Modify transcripts with junctions

We modify the canonical transcripts by introducing the splice junctions.
Then we add the transcript sequence in a fixed-sized window around the
junction positions, the context sequence.

``` r
toy_junc_df
#> # A tibble: 14 × 2
#>    junc_id                    tx_id          
#>    <chr>                      <chr>          
#>  1 chr2_152389996_152392205_- ENST00000409198
#>  2 chr2_152389996_152390729_- ENST00000409198
#>  3 chr2_152389955_152389956_- ENST00000409198
#>  4 chr2_152388410_152392205_- ENST00000409198
#>  5 chr2_152388410_152390729_- ENST00000409198
#>  6 chr2_179415981_179416357_- ENST00000342992
#>  7 chr2_179415987_179415988_- ENST00000342992
#>  8 chr2_179415000_179416357_- ENST00000342992
#>  9 chr2_179445336_179446207_- ENST00000342992
#> 10 chr2_179446225_179446226_- ENST00000342992
#> 11 chr2_179445336_179446633_- ENST00000342992
#> 12 chr2_179642044_179642187_- ENST00000342992
#> 13 chr2_179642146_179642147_- ENST00000342992
#> 14 chr2_179642044_179642431_- ENST00000342992


junc_df <- toy_junc_df %>% 
  add_context_seq(transcripts = toy_transcripts, size = 400, bsg = bsg)


junc_df
#> # A tibble: 14 × 8
#>    junc_id  tx_id  tx_id_alt  junc_pos_tx cts_seq   cts_junc_pos cts_size cts_id
#>    <chr>    <chr>  <chr>            <int> <chr>            <dbl>    <int> <chr> 
#>  1 chr2_15… ENST0… ENST00000…       16412 AAGAAGAC…          200      400 ef606…
#>  2 chr2_15… ENST0… ENST00000…       16517 AAGAAGTA…          200      400 6c189…
#>  3 chr2_15… ENST0… ENST00000…       17290 ACATCTCT…          200      400 c8bd5…
#>  4 chr2_15… ENST0… ENST00000…       16412 AAGAAGAC…          200      400 d41d2…
#>  5 chr2_15… ENST0… ENST00000…       16517 AAGAAGTA…          200      400 db9b3…
#>  6 chr2_17… ENST0… ENST00000…       83789 TGGATTCC…          200      400 744c1…
#>  7 chr2_17… ENST0… ENST00000…       84158 ATTTGAAG…          200      400 5315f…
#>  8 chr2_17… ENST0… ENST00000…       83789 TGGATTCC…          200      400 8eec0…
#>  9 chr2_17… ENST0… ENST00000…       59307 CGGGCTGA…          200      400 5ab65…
#> 10 chr2_17… ENST0… ENST00000…       59288 TTATCTCG…          200      400 c233b…
#> 11 chr2_17… ENST0… ENST00000…       58982 TGGCTATT…          200      400 fddf5…
#> 12 chr2_17… ENST0… ENST00000…        4828 TAGAAGGG…          200      400 ce662…
#> 13 chr2_17… ENST0… ENST00000…        4868 TAGACCTA…          200      400 86af1…
#> 14 chr2_17… ENST0… ENST00000…        4703 GTCTCCTG…          200      400 ec963…
```

### Annotate peptide sequence

Finally, we use the splice junctions to modify the coding sequences
(CDS) of the reference transcripts. The resulting CDS sequences are
translated into protein sequence and further annotated with the peptide
around the junction, the relative position of the splice junction in the
peptide, and the location of the junction in an open reading frame
(ORF).

``` r
junc_df <- junc_df %>% 
  add_peptide(cds=toy_cds, size = 30, bsg = bsg)

junc_df %>% 
  dplyr::select(junc_id, junc_in_orf, peptide_context, peptide_context_junc_pos)
#> # A tibble: 14 × 4
#>    junc_id                    junc_in_orf peptide_context     peptide_context_j…
#>    <chr>                      <lgl>       <chr>                            <dbl>
#>  1 chr2_152389996_152392205_- TRUE        PINRHFKYATQLMNEIC                   14
#>  2 chr2_152389996_152390729_- TRUE        PRHLLAKTAGDQISQIC                   14
#>  3 chr2_152389955_152389956_- FALSE       <NA>                                NA
#>  4 chr2_152388410_152392205_- TRUE        PINRHFKYATQLMNEIKY…                 14
#>  5 chr2_152388410_152390729_- TRUE        PRHLLAKTAGDQISQIKY…                 14
#>  6 chr2_179415981_179416357_- TRUE        PSDPSKFTLAVSPVAGTP…                 14
#>  7 chr2_179415987_179415988_- FALSE       <NA>                                NA
#>  8 chr2_179415000_179416357_- TRUE        PSDPSKFTLAVSPVVPPI…                 14
#>  9 chr2_179445336_179446207_- TRUE        KHYPKDILSKYYQGDST                   14
#> 10 chr2_179446225_179446226_- TRUE        PSDVPDKHYPKDILSKYY…                 14
#> 11 chr2_179445336_179446633_- TRUE        PSDASKAAYARDPQFPPE…                 14
#> 12 chr2_179642044_179642187_- TRUE        TPSDSGEWTVVAQNRLWN…                 14
#> 13 chr2_179642146_179642147_- TRUE        RAGRSSISVILTVEGKMR                  14
#> 14 chr2_179642044_179642431_- TRUE        VVGRPMPETFWFHDAVEH…                 14
```

## Dummy example

In the following a dummy example workflow how to integrate predict
splicing effects from mutations or which detect expressed splice
junctions from RNA-seq data to predict potential neoantigen candidates
with splice2neo. A test case will be added later

``` r
library(splice2neo)
library(tidyverse)
# load genome of choice
library(BSgenome.Hsapiens.UCSC.hg19)
library(AnnotationDbi)

# this is an customized example of a transcript database
# the user can choose the best suited database for their use case
# please find below instruction how to create the database from a gtf file
gtf_url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz"
# parse GTF file as txdb object
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_url)

transcripts <-
  GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)
transcripts_gr <- GenomicFeatures::transcripts(txdb)
# Build a GRangesList with cds composed of individual exon ranges
cds <- GenomicFeatures::cdsBy(txdb, by = c("tx"), use.name = TRUE)

# canonical junctions
# the user can choose the best suited data sets for canonical junctions
# the object `canonical junction` should be a vector of canonical junctions in the junction format
canonical_juncs <-
  c("chr1_33361245_33361511_-",
    "chr1_32849649_32852380_-",
    "chrom_start_end_strand")
# OPTIONAL: if canonical junctions are available as a bed file, the user can transform the bed file into the juntion format:
# for bed files defined by exon-exon boundaries:
canonical_juncs <- bed_to_junc(bed_file = "/path/to/canonical/file.bed", type = "exon-exon")
# for bed files defined by introns:
canonical_juncs <- bed_to_junc(bed_file = "/path/to/canonical/file.bed", type = "intron")


# import RNA data
dat_leafcutter <-
  leafcutter_transform(path = "/your/path/to/leafcutter/results")
# supported events: exon_skip, intron_retention, alt_3prime, alt_5prime, mutex_exons
dat_spladder <-
  spladder_transform(path = "/your/path/to/spladder/results")
dat_rna <-
  generate_combined_dataset(spladder_juncs = dat_spladder, leafcutter_juncs = dat_spladder)

# import & transform SpliceAi results
dat_spliceai <-
  parse_spliceai(vcf_file = "path/to/spliceai/file.vcf")
dat_splicai_formatted <- format_spliceai(dat_spliceai)
dat_spliceai_annotated <-
  annotate_spliceai_junction(var_df = dat_splicai_formatted,
                             transcripts = transcripts,
                             transcripts_gr = transcripts_gr)

# import & transform MMSplice results
dat_mmsplice <- parse_spliceai(infile = "path/to/mmsplice/file.csv")
dat_mmsplice_annotated  <-
  annotate_mmsplice(mmsplice_df = dat_mmsplice, transcripts = transcripts)

# mutation-based junctions
dat_mut <-
  combine_mut_junc(spliceai_juncs = dat_spliceai_annotated, mmsplice_juncs = dat_mmsplice_annotated)

# add information if junction is canonical and if found to be expressed by SplAdder or LeafCutter
dat_mut <- dat_mut %>%
  mutate(is_canonical = is_canonical(junc_id, ref_junc = canonical_juncs, exons_gr = transcripts)) %>%
  mutate(is_in_rnaseq = is_in_rnaseq(junc_id, rna_juncs = dat_rna$junc_id))

# remove canonical junctions for further downstream analysis
dat_for_requantification <- dat_mut %>%
  filter(!is_canonical)

# add context sequences
# a list of GRanges with the transcript needs to be added at the moment
# this will be done within add_context_seq in a future version
dat_for_requantification_cts <- dat_for_requantification %>%
  add_context_seq(size = 400, bsg = BSgenome.Hsapiens.UCSC.hg19, transcripts = tanscripts)


# transform to easyquant-format
dat_easyquant <- dat_for_requantification_cts %>%
  transform_for_requant()
write_delim(dat_easyquant, "path/to/easyquant/input/file.txt", delim = "\t")
# DO RE-QUANTIFICATION WITH EASYQUANT
# https://github.com/TRON-Bioinformatics/easyquant


# add peptide sequence
dat_for_requantification_cts_peptide <-
  dat_for_requantification_cts  %>%
  add_peptide(size = 30, bsg = BSgenome.Hsapiens.UCSC.hg19, cds = cds)

# merge EasyQuant results with data
dat_cts_peptide_requantification <-
  map_requant(path_to_easyquant_folder = "/path/to/easyquant/output_folder",
              junc_tib = dat_for_requantification_cts_peptide)


# EasyQuant results can be imported without direct merging with data
dat_requant <-
  read_requant(path_folder = "/path/to/easyquant/output_folder")
```

## Transcript database

To transform mutations into junction format, a database of transcripts
is required. This database can be created as described below:

``` r
# use gtf file of choice and transform into transcript database
gtf_url <- "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh37_mapping/gencode.v34lift37.annotation.gtf.gz"

# parse GTF file as txdb object
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_url)
saveDb(txdb, file = "/path/to/transripts/txdb.sqlite")

# load 
txdb <- loadDb("/path/to/transripts/txdb.sqlite")
```
