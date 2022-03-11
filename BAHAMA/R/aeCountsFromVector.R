#
#' aeCountsFromVector;
#' Function to convert a vector into a matrix with rows corresponing to the subjects and collums to the PTs
#'
#' @param pts vector with pts ID, rownames are the subject IDs
#'
#' @return A matrix with rows corresponding to the subject, and collumns to the PTs
#'
#' @examples
#' counts = aeCountsFromVector(AE_data)
aeCountsFromVector <- function(pts){
  y_pt <- matrix(0,
                 nrow = length(unique(names(pts))),
                 ncol= length(unique(pts)) )
  rownames(y_pt) = unique(names(pts))
  colnames(y_pt) = unique(pts)
  for(i in 1:length(unique(names(pts)))){
    pts_i <- pts[names(pts) == unique(names(pts))[i]]
    for(j in 1:length(unique(pts_i))){
      y_pt[i, colnames(y_pt) == (unique(pts_i)[j])] <- sum(pts_i == (unique(pts_i)[j]))
    }
  }
  return(y_pt)
}


