#' Generate scale-free network skeleton and simulates corresponding precision matrix
#'
#' A utility function generates scale-free network skeleton and simulates corresponding precision matrix. The non-zero elements of the precision matrix are generated randomly from a uniform distribution with parameters (-upper, -lower) UNION (lower, upper).
#'
#' @importFrom igraph sample_pa
#' @importFrom igraph as_adjacency_matrix
#' @importFrom stats runif
#' @export make_sf_precision_mat
#' @param e Numeric constant, the number of edges to add in each time step, see sample_pa().
#' @param power Numeric constant, the power of the preferential attachment for scale-free network, the default is 1, , see sample_pa().
#' @param p A positive integer, the number of vertices.
#' @param lower A positive number, lower limits of the uniform distribution.
#' @param upper A positive number, upper limits of the uniform distribution.
#' @param diag A small positive number to be added to diagonal elements, which guarantees the precision matrix is positive definite.
#' @returns A precision matrix of size p by p.
#' @examples
#' Omega = make_sf_precision_mat(e=1, p=10)
make_sf_precision_mat=function(e=1, power=1, p=20, lower=.2, upper=.5, diag=0.1){
  g <- sample_pa(n=p, power=power, m=e, directed = FALSE)
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
