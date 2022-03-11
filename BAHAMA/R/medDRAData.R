#' medDRAData object and constructors
#'
#' @slot pt_hlt Mapping between the PTs and HLTs
#' @slot hlt_hlgt Mapping between the HLTs and HLGTs
#' @slot hlgt_soc Mapping between the HLGTs and SOCs
#'
medDRAData <- setClass(Class = "medDRAData",
                       slots = list(
                         pt_hlt = 'matrix',
                         hlt_hlgt = 'matrix',
                         hlgt_soc = 'matrix'
                       )
)


#' medDRAData object and constructors
#'
#' @param medDRAMatrix A matrix of MedDRA IDs. The rows correspond to the MedDRA structure (SOC, HLGT, HLT, PT), each row correspond to a single AE. If multiaxiallity is present, a single AE is represented in multiple rows.
#' @return A medDRAData object
#' @examples
#' medDRADataSet <- medDRADataFromMatrix(tree)
medDRADataFromMatrix <- function(medDRAMatrix){
  if(dim(medDRAMatrix)[1] != 4){
    stop("The rows of the 'medDRAMatrix' have to correspond with the MedDRA levels (SOC, HLGT, HLT, PT)")
  }
  if(length(unique(medDRAMatrix[4,])) < length(unique(medDRAMatrix[3,]))){
    stop("The rows of the 'medDRAMatrix' have to correspond with the MedDRA levels (SOC, HLGT, HLT, PT), place them in the correct order")
  }

  pt_hlt <- matrix(NA, nrow = length(unique(medDRAMatrix[4,])), ncol = 4)
  pt_hlt[,1] <- unique(medDRAMatrix[4,])
  for( i in 1:length(pt_hlt[,1]) ){
    for( j in 1:length(unique(medDRAMatrix[3, medDRAMatrix[4,] == pt_hlt[i,1] ])) )
      pt_hlt[i,j+1] <-  unique(medDRAMatrix[3,medDRAMatrix[4,] == pt_hlt[i,1]])[j]
  }

  hlt_hlgt <- matrix(NA, nrow = length(unique(medDRAMatrix[3,])), ncol = 4)
  hlt_hlgt[,1] <- unique(medDRAMatrix[3,])
  for( i in 1:length(hlt_hlgt[,1]) ){
    for( j in 1:length(unique(medDRAMatrix[2, medDRAMatrix[3,] == hlt_hlgt[i,1] ])) ){
      hlt_hlgt[i,j+1] <-  unique(medDRAMatrix[2, medDRAMatrix[3,] == hlt_hlgt[i,1]])[j]
    }
  }

  hlgt_soc <- matrix(NA, nrow = length(unique(medDRAMatrix[2,])), ncol = 10)
  hlgt_soc[,1] <- unique(medDRAMatrix[2,])
  for( i in 1:length(hlgt_soc[,1]) ){
    for( j in 1:length(unique(medDRAMatrix[1, medDRAMatrix[2,] == hlgt_soc[i,1] ])) ){
      if(j>10) next
      hlgt_soc[i,j+1] <- unique(medDRAMatrix[1, medDRAMatrix[2,] == hlgt_soc[i,1]])[j]
    }
  }
  object <- new("medDRAData", pt_hlt = pt_hlt, hlt_hlgt = hlt_hlgt, hlgt_soc = hlgt_soc)
  return(object)
}


