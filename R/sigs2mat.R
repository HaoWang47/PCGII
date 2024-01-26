#' Utility function for PCGII inference results
#'
#' @description
#' A utility function takes PCGII inference results as input and generates an adjacency matrix corresponding to the significant partial correlations
#'
#' @export sigs2mat
#' @param sigs A dataframe of locations (row, col) of selected edges.
#' @param P A number, the number of nodes in the network.
#' @returns A matrix of size P*(P-1)/2, with 0, 1.
#' @examples
#' edges=cbind.data.frame(row=c(1,2,3,1,6,2,1,6,1,4),
#'                        col=c(2,1,1,3,2,6,6,1,4,1)) # five edges
#' sigs2mat(sigs = edges, P = 6)
sigs2mat=function(sigs, P){
  m=matrix(0,P,P)
  for (h in 1: dim(sigs)[1]){
    m[sigs[h,1],sigs[h,2]]=1
    m[sigs[h,2],sigs[h,1]]=1
  }
  m
}
