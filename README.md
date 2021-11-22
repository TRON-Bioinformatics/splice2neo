
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
#> ✓ ggplot2 3.3.3     ✓ purrr   0.3.4
#> ✓ tibble  3.1.2     ✓ dplyr   1.0.6
#> ✓ tidyr   1.1.3     ✓ stringr 1.4.0
#> ✓ readr   1.4.0     ✓ forcats 0.5.1
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
#> # A tibble: 3 x 1
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

### Modify transcripts with junctions

We modify the canonical transcripts by introducing the splice junctions.
Then we add the transcript sequence in a fixed-sized window around the
junction positions, the context sequence.

``` r
junc_df <- junc_df %>% 
  add_context_seq(size = 400, bsg = bsg)

junc_df %>% 
  select(junc_id, tx_id, junc_pos_tx, cts_seq, cts_junc_pos, cts_id)
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
