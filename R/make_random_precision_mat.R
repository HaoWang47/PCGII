#' Generate unstructured/random network skeleton and simulates corresponding precision matrix
#'
#' @description
#' A utility function generates unstructured/random network skeleton and simulates corresponding precision matrix. The non-zero elements of the precision matrix are generated randomly from a uniform distribution with parameters (-upper, -lower) UNION (lower, upper).
#'
#' @importFrom igraph sample_gnp
#' @importFrom igraph as_adjacency_matrix
#' @importFrom stats runif
#' @export make_random_precision_mat
#' @param eta A number between 0 and 1, the probability for drawing an edge between two arbitrary vertices, i.e. the sparsity of the network.
#' @param p A positive integer, the number of vertices.
#' @param lower A positive number, lower limits of the uniform distribution.
#' @param upper A positive number, upper limits of the uniform distribution.
#' @param diag A small positive number to be added to diagonal elements, which guarantees the precision matrix is positive definite.
#' @returns A precision matrix of size p by p.
#' @examples
#' Omega = make_random_precision_mat(eta=.2, p=10)
make_random_precision_mat=function(eta=.01, p=20, lower=.2, upper=.5, diag=0.1){
  g <- sample_gnp(n=p, p=eta, directed = FALSE)
  omega=as.matrix(as_adjacency_matrix(g))
  for(h1 in 1:(p-1)){
    for(h2 in (h1+1):p){
      if(omega[h1,h2]!=0){
        temp=runif(1, lower, upper)*sample(c(-1,1),size=1)
        omega[h1,h2]=temp
        omega[h2,h1]=temp
      }
    }
  }
  diag(omega)=rowSums(abs(omega)) + diag
  omega
}
