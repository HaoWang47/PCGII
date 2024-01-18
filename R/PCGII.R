#' PCGII() is the function to apply the proposed method to get the estimated partial correlation graph with information incorporation
#'
#' @param df The main expression dataset, an n by p matrix, in which each row corresponds to a sample and each column represents expression/abundance of an omics feature.
#' @param prior The prior set, a k by 2 dataframe, in which each row corresponds to a pair of nodes (any omics features) that are connected under prior belief. Note, prior input has to be dataframe.
#' @param lambda The regularization parameter, used in the node-wise regression. If missing, default lambda will be used which is at the order of sqrt(2*log(p)/n).
#' @returns A list. The list contains estimated partial correlation matrix (Est), sparse partial correlation estimation matrix with threshold (EstThresh), estimated kappa (kappa), estimated test statistics matrix of partial correlations (tscore), sample size (n) and number of nodes (p).
#' @examples
#' # Simulating data
#' library(igraph)
#' library(tidyverse)
#' set.seed(1234567)
#' n=50 # sample size
#' p=30 # number of nodes
#'
#' g=sample_pa(p, power=1, m=1, directed = FALSE) # undirected scale-free network with the power of the preferential attachment set as 1, the number of edges to add in each time step set as 2.
#' plot(g, vertex.size=4, vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.5) # visulize simulated network structure
#' g %>% plot(vertex.size=4, vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.5, layout=layout_in_circle(g))
#' # compute precision matrix structure corresponding to the simulated scale-free networ
#' omega=as_adjacency_matrix(g) %>% as.matrix()
#' for(h1 in 1:(p-1)){
#'   for(h2 in (h1+1):p){
#'     if(omega[h1,h2]!=0){
#'       temp=runif(1, 0.2, 0.5)*sample(c(-1,1),size=1) # randomly assign connection strength, i.e. partial correlations
#'       omega[h1,h2]=temp
#'       omega[h2,h1]=temp
#'     }
#'   }
#' }
#' diag(omega)=rowSums(abs(omega)) # make sure precision matrix is positive definite
#' diag(omega)=diag(omega)+0.10
#'
#' # population covariance matrix, which is used to generate data
#' Sigma=solve(omega)
#' # simulate expression data
#' X = rmvnorm(n = n, sigma = Sigma)
#'
#' lam=2*sqrt(log(p)/n) ## fixed lambda
#'
#' # directed prior network
#' prior_set=matrix(data=c(9,15, 3,4, 5,24, 16,20, 25,22, 28,8, 11,4), nrow=7, ncol=2, byrow = TRUE)
#' colnames(prior_set)=c("row", "col")
#' PCGII_out=PCGII(df=X, prior=as.data.frame(prior_set), lambda = lam)
#' ## Remark: mathematical standardization will be automatically done within the function.
PCGII=function(df, prior, lambda){
  require(glmnet)
  n = dim(df)[1]; p = dim(df)[2]
  t0=2
  IndMatrix = matrix(1, p, p) - diag(rep(1, p))
  Eresidual = matrix(0, n, p) # regression residuals matrix n*p
  CoefMatrix = matrix(0, p, p - 1) # regression coefficient matrix p*p-1
  meanX = colMeans(df)
  X = t(t(df) - meanX)
  XS = matrix(0, n, p)
  # XS: Standardized X
  for (i in 1 : p){
    XS[, i] = X[, i] / sd(X[, i])
  }

  colnames(X)=colnames(df)
  colnames(XS)=colnames(df)

  if(missing(lambda)){
    shat=sqrt(n/(log(p)^3))
    lambda=sqrt(2*(2+0.01)*log(p/shat)/n)
  }

  default_penalty=rep(1,p-1)
  for (i in 1 : p){
    penalty_fac=default_penalty
    temp.node=prior[with(prior,row==i),'col']

    for(nds in temp.node){
      if (nds < i) {penalty_fac[nds]=0} else {penalty_fac[nds-1]=0}
    }

    out = glmnet(XS[, -i], X[, i], family = "gaussian", lambda = lambda, penalty.factor=penalty_fac)
    Coef = out$beta
    CoefMatrix[i, ] = as.numeric(Coef) / apply(X[, -i], 2, sd)

    out = glmnet(XS[, -i], X[, i], family = "gaussian", lambda = lambda)
    Predict = predict(out, XS[, -i], type = "link")
    Eresidual[, i] = X[, i] - Predict
  }

  CovRes = t(Eresidual) %*% Eresidual / n # residuals covariance
  Est = matrix(1, p, p) # estimated partial correlation (rho hat in the paper )

  for (i in 1 : (p - 1)){
    for (j in (i + 1) : p){
      temp = Eresidual[, i] * Eresidual[, j] + Eresidual[, i]^2 * CoefMatrix[j, i] + Eresidual[, j]^2 * CoefMatrix[i, j - 1]
      Est[i, j] = mean(temp) / sqrt(diag(CovRes)[i] * diag(CovRes)[j])
      Est[j, i] = Est[i, j]
    }
  }

  EstThresh = Est * ( abs(Est) >= (t0 * sqrt(log(p) / n) * IndMatrix) )

  kappa = (n / 3) * mean( colSums(Eresidual^4) / (colSums(Eresidual^2))^2 )  # forth moment, a number

  SE=sqrt((kappa*(1-EstThresh^2))^2/n)

  tscore=Est/SE

  return(list(Est=Est,
              tscore=tscore,
              kappa=kappa,
              EstThresh=EstThresh,
              n=n, p=p))

}
