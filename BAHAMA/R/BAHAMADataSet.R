#' BAHAMADataSet object and constructors
#' @slot N_trt The total number of patients in the treatment group
#' @slot N_cont The total number of patients in the control group
#' @slot N_ae The total number of adverse events, after data aggregation
#' @slot N_pt The total number of PTs, after data aggregation
#' @slot N_hlt The total number of HLTs, after data aggregation
#' @slot N_hlt_noPT The total number of HLTs, without a corresponing PT, after data aggregation
#' @slot N_hlt_PT The total number of HLTs, with a corresponing PT, after data aggregation
#' @slot N_hlgt The total number of HLGTs, after data aggregation
#' @slot N_hlgt_noHLT The total number of HLGTs, without a corresponing HLT, after data aggregation
#' @slot N_hlgt_hlt The total number of HLGTs, with a corresponing HLT, after data aggregation
#' @slot N_soc The total number of SOCs, after data aggregation
#' @slot y_hlgt_trt Per HLGTs the incidence in the treatment group, after data aggregation
#' @slot y_hlgt_cont Per HLGTs the incidence in the control group, after data aggregation
#' @slot y_hlt_trt Per HLTs the incidence in the treatment group, after data aggregation
#' @slot y_hlt_cont Per HLTs the incidence in the control group, after data aggregation
#' @slot y_pt_trt Per PTs the incidence in the treatment group, after data aggregation
#' @slot y_pt_cont Per PTs the incidence in the control group, after data aggregation
#' @slot SOC_HLGT The mapping between the HLGTs and SOCs
#' @slot HLGT_HLT The mapping between the HLTs and HLGTs
#' @slot HLT_PT The mapping between the PTs and HLTs
#' @slot N_PT_perHLT The number of PTs per HLT, to correct for data aggregation
#' @slot N_PT_perHLGT The number of PTs per HLGT, to correct for data aggregation
#'
#'
BAHAMADataSet <- setClass(Class = "BAHAMADataSet",
                          slots = list(
                            N_trt = "numeric",
                            N_cont = "numeric",
                            N_ae = "numeric",
                            N_pt = "numeric",
                            N_hlt = "numeric",
                            N_hlt_noPT = "numeric",
                            N_hlt_PT = "numeric",
                            N_hlgt = "numeric",
                            N_hlgt_noHLT = "numeric",
                            N_hlgt_hlt = "numeric",
                            N_soc = "numeric",
                            y_hlt_trt = "numeric",
                            y_hlt_cont = "numeric",
                            y_hlgt_trt= "numeric",
                            y_hlgt_cont = "numeric",
                            y_pt_trt= "numeric",
                            y_pt_cont = "numeric",
                            w_hlgt_soc = "matrix",
                            w_hlt_hlgt = "matrix",
                            w_pt_hlt = "matrix",
                            N_PT_perHLT = "numeric",
                            N_PT_perHLGT = "numeric"
                          )
)

#' BAHAMADataSet object and constructors
#'
#' @param aeCounts A matrix of non-negative integers. Each column of aeCounts correspond to a single AE, each row to a sample.
#' @param BAHAMAData A medDraData object
#' @param sampleData A vector or matrix, with at least 1 column with the treatment group. Each row corresponds to a single sample.
#' @param thresholdIncidence A non-negative interger with the minimal incidence of a single AE
#' @param thresholdStructure A non-negative interger with the minimal number of siblings of a single AE
#' @param alternative_weights A list with a PT weight matrix, a HLT weight matrix and a HLGT weight matrix
#' @return A object to store the input values for the BhaMae method.
#' @examples
#' add(1, 1)
#' add(10, 1)
BAHAMADataSet <- function(aeCounts, medDRAData, sampleData, thresholdIncidence = 5, thresholdStructure = 2, alternative_weights = 0) {
  if(dim(aeCounts)[1] != dim(sampleData)[1]){
    stop("The number of samples in the 'aeCounts' does not equal the number of sample in the 'sampleData'")
  }

  aeCounts <- aeCounts[,order(match(colnames(aeCounts), medDRAData@pt_hlt[,1]))]

  message("Aggregating AEs count based on threshold settings")
  y_pt_trt <- apply(aeCounts * sampleData[,1], 2, sum)
  y_pt_cont <- apply(aeCounts * (1- sampleData[,1]), 2, sum)

  pts_to_include <- aggregating(y_trt = y_pt_trt, y_cont = y_pt_cont,
                                thresholdIncidence = thresholdIncidence, thresholdStructure = thresholdStructure,
                                medDRAData = medDRAData@pt_hlt)

  y_hlt_trt <- c()
  y_hlt_cont <- c()
  for(i in 1:length(pts_to_include$ae_higher_level)){
    y_hlt_trt <- c(y_hlt_trt, sum(y_pt_trt[medDRAData@pt_hlt[,2] == pts_to_include$ae_higher_level[i]]))
    y_hlt_cont <- c(y_hlt_cont, sum(y_pt_cont[medDRAData@pt_hlt[,2] == pts_to_include$ae_higher_level[i]]))
  }
  names(y_hlt_trt) <- pts_to_include$ae_higher_level
  names(y_hlt_cont) <- pts_to_include$ae_higher_level

  hlt_hlgt <- medDRAData@hlt_hlgt[medDRAData@hlt_hlgt[,1] %in% names(y_hlt_trt),]
  hlt_hlgt <- hlt_hlgt[order(match(hlt_hlgt[,1],names(y_hlt_trt))),]

  hlts_to_include <- aggregating(y_trt = y_hlt_trt,
                                 y_cont = y_hlt_cont,
                                 thresholdIncidence = thresholdIncidence,
                                 thresholdStructure = thresholdStructure,
                                 medDRAData = hlt_hlgt)

  y_hlgt_trt <- c()
  y_hlgt_cont <- c()
  for(i in 1:length(hlts_to_include$ae_higher_level)){
    y_hlgt_trt <- c(y_hlgt_trt, sum(y_hlt_trt[hlt_hlgt[,2] == hlts_to_include$ae_higher_level[i]]))
    y_hlgt_cont <- c(y_hlgt_cont, sum(y_hlt_cont[hlt_hlgt[,2] == hlts_to_include$ae_higher_level[i]]))
  }
  names(y_hlgt_trt) <- hlts_to_include$ae_higher_level
  names(y_hlgt_cont) <- hlts_to_include$ae_higher_level

  y_pt_trt <- y_pt_trt[pts_to_include$ae_to_include]
  y_pt_cont <- y_pt_cont[pts_to_include$ae_to_include]
  y_hlt_trt <- y_hlt_trt[hlts_to_include$ae_to_include]
  y_hlt_cont <- y_hlt_cont[hlts_to_include$ae_to_include]

  PT_HLGT <- medDRAData@hlt_hlgt[medDRAData@hlt_hlgt[,1] %in% medDRAData@pt_hlt[pts_to_include$ae_to_include,2],]
  PT_HLGT <- PT_HLGT[order(match(PT_HLGT[,1], unique(medDRAData@pt_hlt[pts_to_include$ae_to_include,2]))),]
  y_hlgt_trt <- y_hlgt_trt[!names(y_hlgt_trt) %in% PT_HLGT[,2]]
  y_hlgt_cont <- y_hlgt_cont[!names(y_hlgt_cont) %in% PT_HLGT[,2]]

  HLGT_HLT <- rbind(PT_HLGT, hlt_hlgt[hlts_to_include$ae_to_include,])

  PT_SOC <- medDRAData@hlgt_soc[medDRAData@hlgt_soc[,1] %in% PT_HLGT[,2],]
  PT_SOC <- PT_SOC[order(match(PT_SOC[,1], unique(PT_HLGT[,2]))),]
  HLT_SOC <- medDRAData@hlgt_soc[medDRAData@hlgt_soc[,1] %in% hlt_hlgt[hlts_to_include$ae_to_include,2],]
  HLT_SOC <- HLT_SOC[order(match(HLT_SOC[,1], unique(hlt_hlgt[hlts_to_include$ae_to_include,]))),]
  SOC_HLGT <- rbind(PT_SOC, HLT_SOC)
  SOC_HLGT <- SOC_HLGT[!duplicated(SOC_HLGT[,1]),]
  SOC_HLGT <- rbind(SOC_HLGT, medDRAData@hlgt_soc[medDRAData@hlgt_soc[,1] %in% names(y_hlgt_trt),])

  N_PT_perHLGT <- medDRAData@hlt_hlgt[ medDRAData@hlt_hlgt[,2] %in% names(y_hlgt_trt),]
  N_PT_perHLGT <- cbind(N_PT_perHLGT, rep(0, dim(N_PT_perHLGT)[1]))
  N_PT_perHLT <- table(medDRAData@pt_hlt[,2:4])
  if(dim(N_PT_perHLGT)[1] != 0){
    for(i in 1:dim(N_PT_perHLGT)[1]){
      N_PT_perHLGT[i,5] <- N_PT_perHLT[names(N_PT_perHLT) == N_PT_perHLGT[i,1]]
    }
  }
  N_PT_perHLT <- N_PT_perHLT[names(N_PT_perHLT) %in% names(hlts_to_include$ae_to_include)[hlts_to_include$ae_to_include] ]
  N_PT_perHLT = as.matrix(N_PT_perHLT)
  N_PT_perHLT = N_PT_perHLT[match(names(y_hlt_trt), rownames(N_PT_perHLT)),]

  if(length(y_hlgt_trt) != 0){
    N_PT_perHLGT_def <- rep(0, length(hlts_to_include$ae_higher_level))
    for(i in 1:length(hlts_to_include$ae_higher_level)){
      N_PT_perHLGT_def[i] <- sum(as.numeric(N_PT_perHLGT[N_PT_perHLGT[,2] == hlts_to_include$ae_higher_level[i],5]))
    }
    names(N_PT_perHLGT_def) <-hlts_to_include$ae_higher_level
    N_PT_perHLGT_def <- N_PT_perHLGT_def[N_PT_perHLGT_def != 0]
    #N_PT_perHLGT_def <- cbind(as.numeric(hlts_to_include$ae_higher_level), N_PT_perHLGT_def)

  }
  else{N_PT_perHLGT_def <- matrix(c(0,0), ncol=2)}

  message("Determine MedDRA weight matrixs")
  if(alternative_weights == 0){
    w_pt_hlt = createWeightMatrixFromMedraMatrix(medDRAData@pt_hlt[pts_to_include$ae_to_include,])
    w_hlt_hlgt = createWeightMatrixFromMedraMatrix(medDRAData@hlt_hlgt)
    w_hlt_hlgt = w_hlt_hlgt[rownames(w_hlt_hlgt) %in% colnames(w_pt_hlt) | rownames(w_hlt_hlgt) %in% names(y_hlt_trt),]
    w_hlt_hlgt = w_hlt_hlgt[c(match(colnames(w_pt_hlt), rownames(w_hlt_hlgt)),
                              match(names(y_hlt_trt), rownames(w_hlt_hlgt))),match(SOC_HLGT[,1], colnames(w_hlt_hlgt))]
    w_hlgt_soc = createWeightMatrixFromMedraMatrix(SOC_HLGT)
    w_hlgt_soc = w_hlgt_soc[match(colnames(w_hlt_hlgt), rownames(w_hlgt_soc)),]
  }
  else{
    if(length(alternative_weights) != 3){
      stop("Alternative weight not correct")
    }
    message("Check alternative weights of PT level")
    w_pt_hlt <- alternative_weights[[1]]
    w_pt_hlt <- w_pt_hlt[pts_to_include$ae_to_include,]
    w_pt_hlt <- w_pt_hlt[, !(colnames(w_pt_hlt) %in% pts_to_include$ae_higher_level)]
    w_pt_hlt <- w_pt_hlt / apply(w_pt_hlt, 1, sum)
    w_pt_hlt <- w_pt_hlt[match(names(y_pt_trt), rownames(w_pt_hlt)),]

    message("Check alternative weights of HLT level")
    w_hlt_hlgt <- alternative_weights[[2]]
    w_hlt_hlgt <- w_hlt_hlgt[rownames(w_hlt_hlgt) %in% colnames(w_pt_hlt) |
                               rownames(w_hlt_hlgt) %in% names(y_hlt_trt), ]
    w_hlt_hlgt <- w_hlt_hlgt[, !(colnames(w_hlt_hlgt) %in%  names(y_hlgt_trt))]
    w_hlt_hlgt <- w_hlt_hlgt / apply(w_hlt_hlgt, 1, sum)
    w_hlt_hlgt <- w_hlt_hlgt[c(match(colnames(w_pt_hlt), rownames(w_hlt_hlgt)),
                               match(names(y_hlt_trt), rownames(w_hlt_hlgt))),]

    message("Check alternative weights of HLGT level")
    w_hlgt_soc <- alternative_weights[[3]]
    w_hlgt_soc <- w_hlgt_soc[rownames(w_hlgt_soc) %in% colnames(w_pt_hlt) |
                               rownames(w_hlgt_soc) %in% names(y_hlt_trt), ]
    w_hlgt_soc <- w_hlgt_soc[, !(colnames(w_hlgt_soc) %in%  names(y_hlgt_trt))]
    w_hlgt_soc <- w_hlgt_soc / apply(w_hlgt_soc, 1, sum)
    w_hlgt_soc <- w_hlgt_soc[c(match(colnames(w_pt_hlt), rownames(w_hlgt_soc)),
                               match(names(y_hlt_trt), rownames(w_hlt_hlgt))),]
  }

  message("Create BAHAMA dataset")
  object <- new("BAHAMADataSet",
                N_trt = sum(sampleData[,1]),
                N_cont = sum(1-sampleData[,1]),
                N_ae = length(y_pt_trt) + length(y_hlt_trt) + length(y_hlgt_trt),
                N_pt = length(y_pt_trt),
                N_hlt = length(unique(medDRAData@pt_hlt[pts_to_include$ae_to_include,2])) + length(y_hlt_trt),
                N_hlt_noPT = length(y_hlt_trt),
                N_hlt_PT = length(unique(medDRAData@pt_hlt[pts_to_include$ae_to_include,2])),
                N_hlgt = dim(w_hlgt_soc)[1],
                N_hlgt_noHLT = length(y_hlgt_trt),
                N_hlgt_hlt = dim(w_hlgt_soc)[1] - length(y_hlgt_trt),
                N_soc = dim(w_hlgt_soc)[2],
                y_hlt_trt = y_hlt_trt,
                y_hlt_cont = y_hlt_cont,
                y_hlgt_trt= y_hlgt_trt,
                y_hlgt_cont = y_hlgt_cont,
                y_pt_trt= y_pt_trt,
                y_pt_cont = y_pt_cont,
                w_hlgt_soc = w_hlgt_soc,
                w_hlt_hlgt = w_hlt_hlgt,
                w_pt_hlt = w_pt_hlt,
                N_PT_perHLT = N_PT_perHLT,
                N_PT_perHLGT = N_PT_perHLGT_def)
  return(object)

}
