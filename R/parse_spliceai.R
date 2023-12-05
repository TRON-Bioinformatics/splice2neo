
#' Parse VCF output file from spliceAI as table
#'
#' requires the vcfR package.
#'
#' @param vcf_file path to a single VCF file, the output of spliceAI
#' @return a [tibble][tibble::tibble-package] with one row per variant
#' annotation. Each input variant can have multiple annotations.
#'
#' @examples
#'
#' spliceai_file <- system.file("extdata", "spliceai_output.vcf", package = "splice2neo")
#' parse_spliceai(spliceai_file)
#'
#' @seealso \code{\link{format_spliceai}}, \code{\link{annotate_mut_effect}}
#' @import dplyr stringr tidyr
#' @export
parse_spliceai <- function(vcf_file){

  # Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL

  col_names <- c("ALLELE","SYMBOL","DS_AG",
                 "DS_AL","DS_DG","DS_DL","DP_AG","DP_AL","DP_DG","DP_DL")

  vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)

  # get fixed fields for variants
  if (nrow(vcf) > 1){
    fix_df <-  vcfR::getFIX(vcf) %>%
      tibble::as_tibble() %>%
      mutate(Key = row_number())
  } else {
    fix_df <-  vcfR::getFIX(vcf) %>%
      tibble::as_tibble_row() %>%
      mutate(Key = row_number())
  }

  # get annotations
  splice_df <- vcf %>%
    vcfR::extract_info_tidy(info_fields = c("SpliceAI")) %>%
    mutate(
      fields = SpliceAI %>% stringr::str_split(",")
    ) %>%
    select(-SpliceAI) %>%
    unnest(fields) %>%
    separate(fields, into = col_names, sep = "\\|") %>%
    mutate(

      # replace missing data with "." with NA
      across(starts_with("DS_"), na_if, "."),
      across(starts_with("DP_"), na_if, "."),

      # convert scores as numeric and positions as integers
      across(starts_with("DS_"), as.numeric),
      across(starts_with("DP_"), as.integer)
    )

  fix_df %>%
    left_join(splice_df, by = "Key")

}



#' Parse VCF output file from SpliceAI with -T flag as table
#'
#' requires the vcfR package.
#'
#' @param vcf_file path to a single VCF file, the output of spliceAI
#' @return a [tibble][tibble::tibble-package] with one row per variant
#' annotation. Each input variant can have multiple annotations.
#'
#' @examples
#'
#' spliceai_file <- system.file("extdata", "spliceai_thresh_output.vcf", package = "splice2neo")
#' parse_spliceai_thresh(spliceai_file)
#'
#' @seealso \code{\link{format_spliceai_thresh}}, \code{\link{annotate_mut_effect}}
#' @import dplyr stringr tidyr
#' @export
parse_spliceai_thresh <- function(vcf_file){

  # read in vcf file
  vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)

  # get fixed fields for variants
  fix_df <-  vcfR::getFIX(vcf) %>%
    tibble::as_tibble() %>%
    mutate(Key = row_number())

  # get annotations
  splice_df <- vcf %>%
    vcfR::extract_info_tidy(info_fields = c("SpliceAI")) %>%
    mutate(fields = SpliceAI %>% stringr::str_split(",")) %>%
    dplyr::select(-SpliceAI) %>%
    tidyr::unnest(fields) %>%
    dplyr::filter(!is.na(fields))
  
  split_rows <- lapply(splice_df$fields, function(x){
    split <- unlist(stringr::str_split(x, "\\|"))
    data.frame(ALLELE = split[1],
               SYMBOL = split[2], 
               splice_sites = paste(split[3:length(split)], collapse = '_'))
  })
  split_rows <- Reduce(rbind, split_rows)
  splice_df  <- cbind(splice_df %>% dplyr::select(-fields), split_rows)
  
  splice_df <- splice_df %>% 
    mutate(fields = splice_sites %>% stringr::str_split("_")) %>%
    dplyr::select(-splice_sites) %>%
    tidyr::unnest(fields) %>%
    tidyr::separate(fields, into = c("effect", "score", "pos_rel"), sep = ":") %>%
    mutate(score = as.numeric(score),
           pos_rel = as.integer(pos_rel), 
           effect = as.factor(effect))

  fix_df %>%
    left_join(splice_df, by = "Key")

}

