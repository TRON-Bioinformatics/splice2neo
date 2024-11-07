
<!-- README.md is generated from README.Rmd. Please edit that file -->

# splice2neo

<!-- badges: start -->

[![R-CMD-check](https://github.com/TRON-Bioinformatics/splice2neo/workflows/R-CMD-check/badge.svg)](https://github.com/TRON-Bioinformatics/splice2neo/actions)
[![Codecov test
coverage](https://codecov.io/gh/TRON-Bioinformatics/splice2neo/branch/master/graph/badge.svg)](https://codecov.io/gh/TRON-Bioinformatics/splice2neo?branch=master)
[![](https://img.shields.io/badge/devel%20version-0.6.12-blue.svg)](https://github.com/TRON-Bioinformatics/splice2neo)
[![](https://img.shields.io/badge/lifecycle-experimental-blue.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![](https://img.shields.io/github/last-commit/TRON-Bioinformatics/splice2neo.svg)](https://github.com/TRON-Bioinformatics/splice2neo/commits/dev)
<!-- badges: end -->

- Documentation: <https://tron-bioinformatics.github.io/splice2neo/>
- Publication: Lang et al. (2024) [Prediction of tumor-specific splicing
  from somatic mutations as a source of neoantigen
  candidates](https://doi.org/10.1093/bioadv/vbae080)

## Overview

This package provides functions for the analysis of splice junctions and
their association with somatic mutations. It integrates the output of
several tools that predict splicing effects from mutations or detect
expressed splice junctions from RNA-seq data into a standardized splice
junction format based on genomic coordinates. Users can filter detected
splice junctions against canonical splice junctions and annotate them
with affected transcript sequences, CDS, and resulting peptide
sequences. Splice2neo supports splice events from alternative 3’/5’
splice sites, exons skipping, intron retentions, exitrons, and mutually
exclusive exons. Integrating splice2neo functions and detection rules
based on splice effect scores and RNA-seq support facilitates the
identification of mutation-associated splice junctions that are
tumor-specific and can encode neoantigen candidates.

Currently, splice2neo supports the output data from the following tools:

- Mutation effect prediction tools:
  - [SpliceAI](https://github.com/Illumina/SpliceAI)
  - [CI-SpliceAI](https://github.com/YStrauch/CI-SpliceAI__Annotation)
  - [MMsplice](https://github.com/gagneurlab/MMSplice_MTSplice)
  - [Pangolin](https://github.com/tkzeng/Pangolin)
- RNA-seq-based splice junction detection tools:
  - [LeafCutter](https://github.com/davidaknowles/leafcutter)
  - [RegTools](https://github.com/griffithlab/regtools)
  - [SplAdder](https://github.com/ratschlab/spladder)
  - [IRfinder](https://github.com/RitchieLabIGH/IRFinder)
  - [SUPPA2](https://github.com/comprna/SUPPA)
  - [STAR](https://github.com/alexdobin/STAR)
  - [StringTie](https://github.com/gpertea/stringtie)

For a more detailed description and a full example workflow see the
[vignette](https://tron-bioinformatics.github.io/splice2neo/articles/splice2neo_workflow.html)

## Installation

The R package splice2neo is not yet on
[CRAN](https://CRAN.R-project.org) or
[Bioconductor](https://www.bioconductor.org/). Therefore, you need to
install splice2neo from
[github](https://github.com/TRON-Bioinformatics/splice2neo).

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("TRON-Bioinformatics/splice2neo")
```

## Example usage

This is a basic example of how to use some functions.

``` r
library(splice2neo)

# load human genome reference sequence
requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
```

### 1 Example data

We start with some example splice junctions provided with the package.

``` r
junc_df <- dplyr::tibble(
  junc_id = toy_junc_id[c(1, 6, 10)]
)

junc_df
#> # A tibble: 3 × 1
#>   junc_id                   
#>   <chr>                     
#> 1 chr2:152389996-152392205:-
#> 2 chr2:179415981-179416357:-
#> 3 chr2:179446225-179446226:-
```

### 2 Add transcripts

Next, we find the transcripts which are in the same genomic region as
the splice junction and that may be affected by the junction.

``` r
junc_df %>% 
  add_tx(toy_transcripts) %>% 
  head()
#> # A tibble: 6 × 3
#>   junc_id                    tx_id           tx_lst      
#>   <chr>                      <chr>           <named list>
#> 1 chr2:152389996-152392205:- ENST00000409198 <GRanges>   
#> 2 chr2:152389996-152392205:- ENST00000172853 <GRanges>   
#> 3 chr2:152389996-152392205:- ENST00000397345 <GRanges>   
#> 4 chr2:152389996-152392205:- ENST00000427231 <GRanges>   
#> 5 chr2:152389996-152392205:- ENST00000618972 <GRanges>   
#> 6 chr2:152389996-152392205:- ENST00000413693 <GRanges>
```

### 3 Modify transcripts with junctions

We modify the canonical transcripts by introducing the splice junctions.
Then we add the transcript sequence in a fixed-sized window around the
junction positions, the context sequence.

``` r
toy_junc_df %>% 
  head()
#> # A tibble: 6 × 2
#>   junc_id                    tx_id          
#>   <chr>                      <chr>          
#> 1 chr2:152389996-152392205:- ENST00000409198
#> 2 chr2:152389996-152390729:- ENST00000409198
#> 3 chr2:152389955-152389956:- ENST00000397345
#> 4 chr2:152388410-152392205:- ENST00000409198
#> 5 chr2:152388410-152390729:- ENST00000409198
#> 6 chr2:179415981-179416357:- ENST00000342992
```

``` r


toy_junc_df %>% 
  add_context_seq(transcripts = toy_transcripts, size = 400, bsg = bsg) %>% 
  head()
#> # A tibble: 6 × 8
#>   junc_id       tx_id tx_mod_id junc_pos_tx cts_seq cts_junc_pos cts_size cts_id
#>   <chr>         <chr> <chr>           <int> <chr>   <chr>           <int> <chr> 
#> 1 chr2:1523899… ENST… ENST0000…       16412 AAGAAG… 200               400 90bfc…
#> 2 chr2:1523899… ENST… ENST0000…       16517 AAGAAG… 200               400 26f77…
#> 3 chr2:1523899… ENST… ENST0000…       21620 AAGAAG… 0,200,1745,…     1945 f1f2c…
#> 4 chr2:1523884… ENST… ENST0000…       16412 AAGAAG… 200               400 d4f9e…
#> 5 chr2:1523884… ENST… ENST0000…       16517 AAGAAG… 200               400 c715a…
#> 6 chr2:1794159… ENST… ENST0000…       83789 TGGATT… 200               400 0128d…
```

### 4 Annotate peptide sequence

Here, we use the splice junctions to modify the coding sequences (CDS)
of the reference transcripts. The resulting CDS sequences are translated
into protein sequence and further annotated with the peptide around the
junction, the relative position of the splice junction in the peptide,
and the location of the junction in an open reading frame (ORF).

``` r

toy_junc_df %>% 
  
  # add peptide sequence
  add_peptide(cds=toy_cds, flanking_size = 13, bsg = bsg) %>% 

  # select subset of columns
  dplyr::select(junc_id, peptide_context, junc_in_orf, cds_description) %>% 
  head()
#> # A tibble: 6 × 4
#>   junc_id                    peptide_context         junc_in_orf cds_description
#>   <chr>                      <chr>                   <lgl>       <chr>          
#> 1 chr2:152389996-152392205:- NRHFKYATQLMNEIC         TRUE        mutated cds    
#> 2 chr2:152389996-152390729:- HLLAKTAGDQISQIC         TRUE        mutated cds    
#> 3 chr2:152389955-152389956:- MLTALYNSHMWSQVMSDGM     TRUE        mutated cds    
#> 4 chr2:152388410-152392205:- NRHFKYATQLMNEIKYRKNYEK… TRUE        mutated cds    
#> 5 chr2:152388410-152390729:- <NA>                    TRUE        truncated cds  
#> 6 chr2:179415981-179416357:- SDPSKFTLAVSPVAGTPDYIDV… TRUE        mutated cds
```

## Issues

Please report issues here:
<https://github.com/TRON-Bioinformatics/splice2neo/issues>

## Reference

Lang F, Sorn P, Suchan M, Henrich A, Albrecht C, Köhl N, Beicht A,
Riesgo-Ferreiro P, Holtsträter C, Schrörs B, Weber D, Löwer M, Sahin U,
Ibn-Salem J. Prediction of tumor-specific splicing from somatic
mutations as a source of neoantigen candidates. Bioinform Adv. 2024 May
29;4(1):vbae080. doi:
[10.1093/bioadv/vbae080](https://doi.org/10.1093/bioadv/vbae080). PMID:
38863673; PMCID: PMC11165244.
