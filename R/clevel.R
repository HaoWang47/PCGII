#' Get the estimated partial correlation graph without information incorporation
#'
#' @description
#' clevel() is the function to apply the method originally proposed in paper "Qiu, Y., & Zhou, X. H. (2020). Estimating c-level partial correlation graphs with application to brain imaging". It is used to get the estimated partial correlation graph without information incorporation. Remark: mathematical standardization will be automatically done within the function.
#'
#' @importFrom stats sd
#' @importFrom glmnet glmnet
#' @importFrom stats predict
#' @export clevel
#' @param df The main expression dataset, an n by p matrix, in which each row corresponds to a sample and each column represents expression/abundance of an omics feature.
#' @param lambda The regularization parameter, used in the node-wise regression. If missing, default lambda will be used which is at the order of sqrt(2*log(p)/n).
#' @returns A list. The list contains estimated partial correlation matrix (Est), sparse partial correlation estimation matrix with threshold (EstThresh), estimated kappa (kappa), estimated test statistics matrix of partial correlations (tscore), sample size (n) and number of nodes (p).
#' @examples
#' library(PCGII)
#' library(corpcor)
#' library(glmnet)
#' library(igraph)
#' library(Matrix)
#' library(mvtnorm)
#' # Simulating data
#' set.seed(1234567)
#' n=50 # sample size
#' p=30 # number of nodes
#'
#' Omega=make_random_precision_mat(eta=.01, p=p)
#'
#' # population covariance matrix, which is used to generate data
#' Sigma=solve(Omega)
#' # simulate expression data
#' X = rmvnorm(n = n, sigma = Sigma)
#'
#' lam=2*sqrt(log(p)/n) ## fixed lambda
#'
#' CLEVEL_out=clevel(df=X, lambda = lam)
#' inference_out=inference(list=CLEVEL_out)
#' diag(inference_out)=0
#' net=graph_from_adjacency_matrix(inference_out, mode = "undirected")
#'    plot(net,
#'         vertex.size=4,
#'         vertex.label.dist=0.5,
#'         vertex.color="red",
#'         edge.arrow.size=0.5,
#'         layout=layout_in_circle(net))
clevel=function(df, lambda){
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

  for (i in 1 : p){
    out = glmnet(XS[, -i], X[, i], family = "gaussian", lambda = lambda)

    Coef = out$beta
    Predict = predict(out, XS[, -i], type = "link")
    CoefMatrix[i, ] = as.numeric(Coef) / apply(X[, -i], 2, sd)
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

  # sparse partial correlation estimation with threshold (rho ~ in the paper)
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
