#' aggregating
#'
#' @param y_trt
#' @param y_cont
#' @param thresholdIncidence
#' @param thresholdStructure
#' @param medDRAData
#'
aggregating <- function(y_trt, y_cont, thresholdIncidence, thresholdStructure, medDRAData){
  ae_threshold_incidence <-  (y_trt + y_cont)>(thresholdIncidence-1)
  freq_per_higher_level <- table(medDRAData[ae_threshold_incidence,c(2:(dim(medDRAData)[2]))])
  ae_under_threshold_structure <- names(freq_per_higher_level)[freq_per_higher_level<thresholdStructure]
  ae_threshold_strucutre_pri <- !medDRAData[,2] %in% ae_under_threshold_structure
  ae_threshold_strucutre_sec <- !medDRAData[,3] %in% ae_under_threshold_structure

  ae_to_include <- ae_threshold_incidence & ae_threshold_strucutre_pri & ae_threshold_strucutre_sec

  ae_higher_level <- medDRAData[(medDRAData[,2] %in% ae_under_threshold_structure),]
  ae_excluded <- medDRAData[!ae_threshold_incidence,]
  if(is.null(dim(ae_excluded))){
    ae_excluded<- ae_excluded[!ae_excluded[2] %in% medDRAData[ae_threshold_incidence,2]]
  }
  else {ae_excluded<- ae_excluded[!ae_excluded[,2] %in% medDRAData[ae_threshold_incidence,2],]}
  if(is.null(dim(ae_excluded))){
    ae_higher_level <- unique(c(ae_higher_level[,2] , ae_excluded[2]))
  }
  else {ae_higher_level <- unique(c(ae_higher_level[,2] , ae_excluded[,2]))}


  return(list(ae_to_include = ae_to_include,
              ae_higher_level = ae_higher_level))
}
