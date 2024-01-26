#' Generate block-diagonal matrix of size p by p
#'
#' @description
#' A utility function generates block-diagonal matrix of size p by p with blocks B1, B2, ..., Bk. Each block matrix is of size blocksize by blocksize. The off-diagonal elements in block matrix are generated from uniform (min.beta, max.beta). The diagonal elements in block matrix are generated from uniform (1, 1.25).
#'
#' @importFrom Matrix bdiag
#' @importFrom corpcor is.positive.definite
#' @importFrom stats runif
#' @export makeBlockDiag
#' @param blocksize A positive integer, the dimension of the block matrix. Note, 'blocksize' has to be a factor of 'p'.
#' @param p A positive integer, the size of the block-diagonal matrix.
#' @param min.beta A positive number, lower limits of the uniform distribution.
#' @param max.beta A positive number, upper limits of the uniform distribution.
#' @returns A block-diagonal matrix of size 'p' by 'p'.
#' @examples
#' mat = makeBlockDiag(blocksize=4, p=20)
makeBlockDiag=function(blocksize=4, p=20, min.beta=0.3, max.beta=0.9){
  # blocksize has to be a factor of p
  reps=p/blocksize
  S=list()
  for (i in 1:reps) {
    bd=matrix(runif(1, min.beta, max.beta), blocksize, blocksize)
    diag(bd)=runif(1,1,1.25)
    while(!is.positive.definite(bd)){diag(bd)=diag(bd)+0.01}
    S[[i]]=bd
  }
  as.matrix(Matrix::bdiag(S))
}
