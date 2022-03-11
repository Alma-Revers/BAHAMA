#' Create weight matrix from MedDRA matrix
#'
#' @param medDRAmapping
#'
#' @return A weight matrix, every row sums up to 1
#'
#' @examples
#' w = createWeightMatrixFromMedraMatrix(tree[1:2,])
createWeightMatrixFromMedraMatrix <- function(medDRAmapping){
  max_N <- dim(medDRAmapping)[2]
  parent = unique(c(medDRAmapping[,2:max_N]))
  parent = parent[!is.na(parent)]
  w = matrix(0, nrow = dim(medDRAmapping)[1], ncol = length(parent),
             dimnames = list(medDRAmapping[,1],
                             parent))
  for(i in 1:dim(medDRAmapping)[1]){
    parents_i = medDRAmapping[i,2:max_N]
    parents_i = parents_i[!is.na(parents_i)]
    w[i, colnames(w) %in% parents_i] = 1/length(parents_i)
  }
  return(w)
}
