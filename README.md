
<!-- README.md is generated from README.Rmd. Please edit that file -->

# splice2neo

<!-- badges: start -->

[![R-CMD-check](https://github.com/TRON-Bioinformatics/splice2neo/workflows/R-CMD-check/badge.svg)](https://github.com/TRON-Bioinformatics/splice2neo/actions)
[![Codecov test
coverage](https://codecov.io/gh/TRON-Bioinformatics/splice2neo/branch/master/graph/badge.svg)](https://codecov.io/gh/TRON-Bioinformatics/splice2neo?branch=master)
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
require(tidyverse)
#> Loading required package: tidyverse
#> ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──
#> ✓ ggplot2 3.3.5     ✓ purrr   0.3.4
#> ✓ tibble  3.1.6     ✓ dplyr   1.0.7
#> ✓ tidyr   1.1.4     ✓ stringr 1.4.0
#> ✓ readr   2.0.1     ✓ forcats 0.5.1
#> ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
```

### Example data

We start with some example splice junctions provided with the package.

``` r
junc_df <- tibble(
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
junc_df <- junc_df %>% 
  add_context_seq(size = 400, bsg = bsg)

junc_df %>% 
  select(junc_id, tx_id, junc_pos_tx, cts_seq, cts_junc_pos, cts_id)
#> # A tibble: 21 × 6
#>    junc_id                    tx_id  junc_pos_tx cts_seq    cts_junc_pos cts_id 
#>    <chr>                      <chr>        <int> <chr>             <dbl> <chr>  
#>  1 chr2_152389996_152392205_- ENST0…       16412 AAGAAGACT…          199 ef6060…
#>  2 chr2_152389996_152392205_- ENST0…       16412 AAGAAGACT…          199 ef6060…
#>  3 chr2_152389996_152392205_- ENST0…       21515 AAGAAGACT…          199 729100…
#>  4 chr2_152389996_152392205_- ENST0…       21515 AAGAAGACT…          199 ef6060…
#>  5 chr2_152389996_152392205_- ENST0…       21515 AAGAAGACT…          199 729100…
#>  6 chr2_152389996_152392205_- ENST0…        5502 AAGAAGACT…          199 ef6060…
#>  7 chr2_152389996_152392205_- ENST0…       21312 AAGAAGACT…          199 729100…
#>  8 chr2_152389996_152392205_- ENST0…       21312 AAGAAGACT…          199 ef6060…
#>  9 chr2_152389996_152392205_- ENST0…         576 AAGAAGACT…          199 8c2b82…
#> 10 chr2_179415981_179416357_- ENST0…       83789 TGGATTCCA…          199 744c11…
#> # … with 11 more rows
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
  mutate(cds_lst = as.list(toy_cds[tx_id])) %>% 
  add_peptide(size = 30, bsg = bsg)

junc_df %>% 
  select(junc_id, tx_id, junc_in_orf, peptide_context, peptide_context_junc_pos)
#> # A tibble: 21 × 5
#>    junc_id                    tx_id junc_in_orf peptide_context peptide_context…
#>    <chr>                      <chr> <lgl>       <chr>                      <dbl>
#>  1 chr2_152389996_152392205_- ENST… TRUE        PINRHFKYATQLMN…               14
#>  2 chr2_152389996_152392205_- ENST… TRUE        PINRHFKYATQLMN…               14
#>  3 chr2_152389996_152392205_- ENST… TRUE        PINRHFKYATQLMN…               14
#>  4 chr2_152389996_152392205_- ENST… TRUE        PINRHFKYATQLMN…               14
#>  5 chr2_152389996_152392205_- ENST… TRUE        PINRHFKYATQLMN…               14
#>  6 chr2_152389996_152392205_- ENST… TRUE        PINRHFKYATQLMN…               14
#>  7 chr2_152389996_152392205_- ENST… TRUE        PINRHFKYATQLMN…               14
#>  8 chr2_152389996_152392205_- ENST… TRUE        PINRHFKYATQLMN…               14
#>  9 chr2_152389996_152392205_- ENST… TRUE        PINRHFKYATQLMN…               14
#> 10 chr2_179415981_179416357_- ENST… TRUE        PSDPSKFTLAVSPV…               14
#> # … with 11 more rows
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
txdb <- loadDb("/path/to/transripts/txdb.sqlite")
transcripts <-
  GenomicFeatures::exonsBy(txdb, by = c("tx"), use.names = TRUE)
transcripts_gr <- GenomicFeatures::transcripts(txdb)
# Build a GRangesList with cds composed of individual exon ranges
cds <- GenomicFeatures::cdsBy(txdb, by = c("tx"), use.name = TRUE)

# canonical junctions
# the user can choose the best suited datasets for canoical junctions
# the object `canonical junction` should be a vector of canonical junctions in the junction format
canonical_juncs <-
  c("chr1_33361245_33361511_-",
    "chr1_32849649_32852380_-",
    "chrom_start_end_strand")


# import RNA data
dat_leafcutter <-
  leafcutter_transform(path = "/your/path/to/leafcutter/results")
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
  mutate(tx_lst = as.list(transcripts[tx_id])) %>%
  add_context_seq(size = 400, bsg = BSgenome.Hsapiens.UCSC.hg19)


# transform to easyquant-format
dat_easyquant <- dat_for_requantification_cts %>%
  transform_for_requant()
write_delim(dat_easyquant, "path/to/easyquant/input/file.txt", delim = "\t")
# DO RE-QUANTIFICATION WITH EASYQUANT


# add peptide sequence
# a list of GRanges with the CDS needs to be added at the moment
# this will be done within add_peptide in a future version
dat_for_requantification_cts_peptide <-
  dat_for_requantification_cts  %>%
  mutate(cds_lst = as.list(cds[tx_id])) %>%
  add_peptide(size = 30, bsg = BSgenome.Hsapiens.UCSC.hg19)

# merge EasyQuant results with data
dat_cts_peptide_requantification <-
  map_requant(path_to_easyquant_folder = "/path/to/easyuant/output_folder",
              junc_tib = dat_for_requantification_cts_peptide)


# EasyQuant results can be imported without direct merging with data
dat_requant <-
  read_requant(path_folder = "/path/to/easyuant/output_folder")
```

## Transcript database

To transform mutations into junction format, a database of transcripts
is required. This database can be created as described below:

``` r
# use gtf file of choice and transform into transcript database
gtf_file = "/path/to/human/gencode/gencode.annotation.gtf"
# parse GTF file as txdb object
txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file)
saveDb(txdb, file = "/path/to/transripts/txdb.sqlite")
```
