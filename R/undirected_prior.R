#' Pre-process the input prior set to ensure the input prior set corresponds to an undirected prior network
#'
#' @description
#' An utility function to pre-process the input prior set. This function will ensure the input prior set corresponds to an undirected prior network. If the prior network is believed to be directed, no pre-processing of the prior set is needed. Remark: this function is not necessary. Prior set should be considered carefully before running the network analysis. If the prior network connections are believed to be undirected while the prior set only includes one way connections for simplicity, this function will duplicate the connections and swap the direction automatically.
#'
#' @importFrom dplyr arrange
#' @export undirected_prior
#' @param prior A k by 2 data.frame of prior set, in which each row corresponds to a pair of nodes (any omics features) that are connected under prior belief
#' @returns A 2-column data.frame of pre-processed prior set, in which the connection between any pair of nodes is undirected.
#' @examples
#' prior=as.data.frame(matrix(c(1,2,1,5,1,10), ncol=2, byrow=TRUE))
#' ## a prior set of 3 connections (1-2, 1-3, 1-10)
#' colnames(prior)=c("row", "col")
#' undirected=undirected_prior(prior)
undirected_prior = function(prior){
  unique(arrange(rbind(prior, transform(prior, row = pmax(row, col), col = pmin(row, col))), row, col))
}
