
#' Combines tibbles with junctions from SpliceAI and MMsplice
#'
#' @param spliceai_juncs A tibble the junctions identified by SpliceAI
#' @param mmsplice_juncs A tibble the junctions identified by MMsplice in
#'   standardized format
#'
#' @return A combined table with unique junctions. The columns DNA_tool
#'   contains information which tools identified the given junction
#'
#'
#' @import dplyr
#' @export
combine_mut_junc <- function(spliceai_juncs, mmsplice_juncs){

  spliceai_juncs <- spliceai_juncs %>%
    mutate(junc_tx_id = paste0(junc_id, "_", tx_id))%>%
    dplyr::rename(score_spliceai = score)%>%
    dplyr::rename(class_spliceai = class)

  mmsplice_juncs <- mmsplice_juncs %>%
    mutate(junc_tx_id = paste0(junc_id, "_", transcript_id),
           tx_id = transcript_id)%>%
    dplyr::rename(delta_logit_psi_mmsplice = delta_logit_psi) %>%
    dplyr::rename(pathogenicity_mmsplice = pathogenicity)%>%
    dplyr::rename(class_mmsplice = effect)%>%
    dplyr::rename(efficiency_mmsplice = efficiency)

  #juncs <- bind_rows(spliceai_juncs, mmsplice_juncs)
  juncs <- spliceai_juncs %>%
    full_join(mmsplice_juncs, by = "junc_tx_id" )
  juncs <- juncs %>%
    distinct(junc_tx_id, .keep_all = T) %>%
    mutate(identified_by_spliceai = ifelse(junc_tx_id %in% spliceai_juncs$junc_tx_id, TRUE, FALSE),
           identified_by_mmsplice = ifelse(junc_tx_id %in% mmsplice_juncs$junc_tx_id, TRUE, FALSE )) %>%
    mutate(junc_id = ifelse(!is.na(junc_id.x), junc_id.x, junc_id.y))%>%
    mutate(mut_id = ifelse(!is.na(mut_id.x), mut_id.x, mut_id.y))%>%
    mutate(tx_id = ifelse(!is.na(tx_id.x), tx_id.x, tx_id.y))%>%
    dplyr::select(junc_tx_id, mut_id, tx_id, junc_id, junc_tx_id,
           identified_by_spliceai, identified_by_mmsplice,
           score_spliceai, class_spliceai, delta_logit_psi_mmsplice,
           pathogenicity_mmsplice, efficiency_mmsplice,  class_mmsplice)
  return(juncs)

}
