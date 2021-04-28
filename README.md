
<!-- README.md is generated from README.Rmd. Please edit that file -->

# splice2neo

<!-- badges: start -->

<!-- badges: end -->

The goal of splice2neo is to …

## Installation

### Install this package from this GitLab

This R package is not yet on [CRAN](https://CRAN.R-project.org) or
[Bioconductor](https://www.bioconductor.org/). Therefore, you have to
install it form this GitLab repository.

However, this repository is a private GitLab reposoritory and therefore
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
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```
