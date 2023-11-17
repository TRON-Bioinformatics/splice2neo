
#' Parse VCF output file from CI-SpliceAI with -t flag as table
#'
#' requires the vcfR package.
#'
#' @param vcf_file path to a single VCF file, the output of spliceAI
#' @return a [tibble][tibble::tibble-package] with one row per variant
#' annotation. Each input variant can have multiple annotations.
#'
#' @examples
#'
#' cispliceai_file <- system.file("extdata", "cispliceai_thresh_output.vcf", package = "splice2neo")
#' parse_cispliceai(cispliceai_file)
#'
#' @seealso \code{\link{format_cispliceai}}, \code{\link{annotate_mut_effect}}
#' @import dplyr stringr tidyr
#' @export
parse_cispliceai_thresh <- function(vcf_file){
  
  # read in vcf file
  vcf <- vcfR::read.vcfR(vcf_file, verbose = FALSE)
  
  # get fixed fields for variants
  fix_df <-  vcfR::getFIX(vcf) %>%
    tibble::as_tibble() %>%
    mutate(Key = row_number())
  
  # get annotations
  splice_df <- vcf %>%
    vcfR::extract_info_tidy(info_fields = c("CI-SpliceAI")) %>%
    mutate(fields = `CI-SpliceAI` %>% stringr::str_split(",")) %>%
    dplyr::select(-`CI-SpliceAI`) %>%
    tidyr::unnest(fields) %>%
    dplyr::filter(!is.na(fields))
  
  split_rows <- lapply(splice_df$fields, function(x){
    split <- unlist(stringr::str_split(x, "\\|"))
    if(length(split) == 2){
      data.frame(ALLELE = split[1],
                 SYMBOL = split[2], 
                 splice_sites = NA)
    } else {
      data.frame(ALLELE = split[1],
                 SYMBOL = split[2], 
                 splice_sites = paste(split[3:length(split)], collapse = '_'))
    }
  })
  split_rows <- Reduce(rbind, split_rows)
  splice_df  <- cbind(splice_df %>% dplyr::select(-fields), split_rows)
  
  splice_df <- splice_df %>% 
    dplyr::filter(!is.na(splice_sites)) %>%
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
