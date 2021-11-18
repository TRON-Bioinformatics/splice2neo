
<!-- README.md is generated from README.Rmd. Please edit that file -->

# splice2neo

<!-- badges: start -->

[![R-CMD-check](https://github.com/TRON-Bioinformatics/splice2neo/workflows/R-CMD-check/badge.svg)](https://github.com/TRON-Bioinformatics/splice2neo/actions)
<!-- badges: end -->

This package provides functions for the analysis of alternative or
aberrant splicing junctions and their creation from or association with
somatic mutations. It integrates the output of several tools which
predict splicing effects from mutation (SpliceAI, MMSplice) or which
detect expressed splice junctions from RNA-seq data (Leafcutter,
Spladder) into a common splice junction format based on genomic
coordinates. Splice junctions can be annotated with affected transcript
sequences, CDS, and resulting peptide sequences. Resulting peptide
sequences may form neoantigens.

Website: <https://tron-bioinformatics.github.io/splice2neo/>

## Installation

### Install this package from this GitLab

This R package is not yet on [CRAN](https://CRAN.R-project.org) or
[Bioconductor](https://www.bioconductor.org/). Therefore, you have to
install it form this GitHub repository.

``` r
install.packages("devtools") # if needed

devtools::install_github("TRON-Bioinformatics/splice2neo")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(splice2neo)

# load human genome reference sequence
requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
```

### Start with some example splice junctions

``` r
# load toy example splice junctions
junc_df <- tibble::tibble(
  junc_id = toy_junc_id[c(1, 6, 10)]
)

junc_df
#> # A tibble: 3 x 1
#>   junc_id                   
#>   <chr>                     
#> 1 chr2_152389996_152392205_-
#> 2 chr2_179415981_179416357_-
#> 3 chr2_179446225_179446226_-
```

### Add transcripts that are affected by the splice junction

``` r
junc_df <- junc_df %>% 
  add_tx(toy_transcripts)

junc_df
#> # A tibble: 21 x 3
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

### Modify canonical transcripts by altered splice junction and add context sequence around junction position

``` r
junc_df <- junc_df %>% 
  add_context_seq(size = 400, bsg = bsg)

junc_df %>% 
  dplyr::select(junc_id, tx_id, junc_pos_tx, cts_seq, cts_junc_pos, cts_id)
#> # A tibble: 21 x 6
#>    junc_id     tx_id   junc_pos_tx cts_seq               cts_junc_pos cts_id    
#>    <chr>       <chr>         <int> <chr>                        <dbl> <chr>     
#>  1 chr2_15238… ENST00…       16412 AAGAAGACTTGACTTGGCTT…          199 ef6060403…
#>  2 chr2_15238… ENST00…       16412 AAGAAGACTTGACTTGGCTT…          199 ef6060403…
#>  3 chr2_15238… ENST00…       21515 AAGAAGACTTGACTTGGCTT…          199 729100c15…
#>  4 chr2_15238… ENST00…       21515 AAGAAGACTTGACTTGGCTT…          199 ef6060403…
#>  5 chr2_15238… ENST00…       21515 AAGAAGACTTGACTTGGCTT…          199 729100c15…
#>  6 chr2_15238… ENST00…        5502 AAGAAGACTTGACTTGGCTT…          199 ef6060403…
#>  7 chr2_15238… ENST00…       21312 AAGAAGACTTGACTTGGCTT…          199 729100c15…
#>  8 chr2_15238… ENST00…       21312 AAGAAGACTTGACTTGGCTT…          199 ef6060403…
#>  9 chr2_15238… ENST00…         576 AAGAAGACTTGACTTGGCTT…          199 8c2b828f5…
#> 10 chr2_17941… ENST00…       83789 TGGATTCCATGTTGAAAAGA…          199 744c11d66…
#> # … with 11 more rows
```

### Add resulting CDS and peptide sequence

``` r
junc_df <- junc_df %>% 
  dplyr::mutate(cds_lst = as.list(toy_cds[tx_id])) %>% 
  add_peptide(size = 30, bsg = bsg)

junc_df %>% 
  dplyr::select(junc_id, tx_id, junc_in_orf, peptide_context, peptide_context_junc_pos)
#> # A tibble: 21 x 5
#>    junc_id        tx_id     junc_in_orf peptide_context       peptide_context_j…
#>    <chr>          <chr>     <lgl>       <chr>                              <dbl>
#>  1 chr2_15238999… ENST0000… TRUE        PINRHFKYATQLMNEIC                     14
#>  2 chr2_15238999… ENST0000… TRUE        PINRHFKYATQLMNEIC                     14
#>  3 chr2_15238999… ENST0000… TRUE        PINRHFKYATQLMNEIC                     14
#>  4 chr2_15238999… ENST0000… TRUE        PINRHFKYATQLMNEIC                     14
#>  5 chr2_15238999… ENST0000… TRUE        PINRHFKYATQLMNEIC                     14
#>  6 chr2_15238999… ENST0000… TRUE        PINRHFKYATQLMNEIC                     14
#>  7 chr2_15238999… ENST0000… TRUE        PINRHFKYATQLMNEIC                     14
#>  8 chr2_15238999… ENST0000… TRUE        PINRHFKYATQLMNEIC                     14
#>  9 chr2_15238999… ENST0000… TRUE        PINRHFKYATQLMNEIC                     14
#> 10 chr2_17941598… ENST0000… TRUE        PSDPSKFTLAVSPVAGTPDY…                 14
#> # … with 11 more rows
```
