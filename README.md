
<!-- README.md is generated from README.Rmd. Please edit that file -->

# splice2neo

<!-- badges: start -->
<!-- badges: end -->

This package provides functions for the analysis of alternative or
aberrant splicing junctions and their creation from or association with
somatic mutations. It integrates the output of several tools which
predict splicing effects from mutation or RNA-seq data into a common
splice junction format based on genomic coordinates. Splice junctions
can be annotated with affected transcript sequences, CDS, and resulting
peptide sequences.

Website: <https://tron.pages.gitlab.rlp.net/splice2neo>

## Installation

### Install this package from this GitLab

This R package is not yet on [CRAN](https://CRAN.R-project.org) or
[Bioconductor](https://www.bioconductor.org/). Therefore, you have to
install it form this GitLab repository.

However, this repository is a private GitLab repository and therefore
you have to create an personal access token (PAT) first. This is
described
[here](https://docs.gitlab.com/ee/user/profile/personal_access_tokens.html).
As Scope use `api` or `read_api`.

``` r
install.packages("remotes") # if needed

remotes::install_gitlab("tron/splice2neo", host = "gitlab.rlp.net", auth_token = "YOUR_ACCESS_TOKEN")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(splice2neo)

# load human genome reference sequence
requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE)
bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

# load toy example splice junctions
junc_df <- tibble::tibble(
  junc_id = toy_junc_id[c(1, 6, 10)]
)

junc_df
#> # A tibble: 3 × 1
#>   junc_id                   
#>   <chr>                     
#> 1 chr2_152389996_152392205_-
#> 2 chr2_179415981_179416357_-
#> 3 chr2_179446225_179446226_-
# add transcripts that are affected by splice junction
junc_df <- junc_df %>% 
  add_tx(toy_transcripts)

# modify transcripts by junctions and add context sequence around junction position
junc_df <- junc_df %>% 
  add_context_seq(size = 400, bsg = bsg)

junc_df
#> # A tibble: 21 × 10
#>    junc_id  tx_id  tx_lst tx_alt_lst tx_id_alt junc_pos_tx cts_seq  cts_junc_pos
#>    <chr>    <chr>  <name> <named li> <chr>           <int> <chr>           <dbl>
#>  1 chr2_15… ENST0… <GRan… <GRanges>  ENST0000…       16412 AAGAAGA…          199
#>  2 chr2_15… ENST0… <GRan… <GRanges>  ENST0000…       16412 AAGAAGA…          199
#>  3 chr2_15… ENST0… <GRan… <GRanges>  ENST0000…       21515 AAGAAGA…          199
#>  4 chr2_15… ENST0… <GRan… <GRanges>  ENST0000…       21515 AAGAAGA…          199
#>  5 chr2_15… ENST0… <GRan… <GRanges>  ENST0000…       21515 AAGAAGA…          199
#>  6 chr2_15… ENST0… <GRan… <GRanges>  ENST0000…        5502 AAGAAGA…          199
#>  7 chr2_15… ENST0… <GRan… <GRanges>  ENST0000…       21312 AAGAAGA…          199
#>  8 chr2_15… ENST0… <GRan… <GRanges>  ENST0000…       21312 AAGAAGA…          199
#>  9 chr2_15… ENST0… <GRan… <GRanges>  ENST0000…         576 AAGAAGA…          199
#> 10 chr2_17… ENST0… <GRan… <GRanges>  ENST0000…       83789 TGGATTC…          199
#> # … with 11 more rows, and 2 more variables: cts_size <int>, cts_id <chr>
```
