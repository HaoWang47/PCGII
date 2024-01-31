#' Conduct simultaneous inference of estimated partial correlations
#'
#' @description
#' Inference() is the function to conduct simultaneous inference of estimated partial correlations.
#'
#' @export inference
#' @importFrom stats pnorm
#' @param list A list returned by either `PCGII()` or `clevel()`.
#' @param alpha A pre-determined False Discovery Rate. Nominal FDR is set at 0.05 by default.
#' @returns An adjacency matrix of significant partial correlations.
#' @examples
#' library(igraph)
#' library(PCGII)
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
#' # directed prior network
#' prior_set=as.data.frame(matrix(data=c(5,6, 28,24), nrow=2, ncol=2, byrow = TRUE))
#' colnames(prior_set)=c("row", "col")
#' prior_set=undirected_prior(prior_set)
#' PCGII_out=PCGII(df=X, prior=prior_set, lambda = lam)
#' inference_out=inference(list=PCGII_out)
#' diag(inference_out)=0
#' net=graph_from_adjacency_matrix(inference_out, mode = "undirected")
#'    plot(net, vertex.size=4,
#'         vertex.label.dist=0.5,
#'         vertex.color="red",
#'         edge.arrow.size=0.5,
#'         layout=layout_in_circle(net))
inference=function(list, alpha=0.05){
  Est=list$Est
  tscore=list$tscore
  kappa=list$kappa
  EstThresh=list$EstThresh
  n=list$n; p=list$p;
  t0=2; tau = seq(0, 3.5, 0.01); smax = n / 2; lentau = length(tau)


  resprop = list() # selected edges with different tau's, a list of 351 elements
  rejectprop = c()
  for (i in 1 : lentau){ # tau vary from 0 to 3.50 by 0.01, length=351
    Threshold = tau[i] * sqrt(kappa * log(p) / n) * (1 - EstThresh^2)

    # c=0
    SRec = 1 * (abs(Est) > Threshold) # selected edge (matrix with 0 & 1) at tau[i]
    NoNSRec = 1 * (SRec == 0)
    resprop[[i]] = which(SRec == 1, arr.ind = TRUE) # select those significant edges at tau[i], off-diagonal elements, first columns, then second columns
    rejectprop = c(rejectprop, max(1, (sum(SRec) - p)))
  }

  # c=0
  FDPprop = 2 * (p * (p - 1)) * ( 1 - pnorm( tau * sqrt(log(p)) ) ) / rejectprop # FDP corresponding to each tau (page 10)

  FDPresprop = c()

  # determine thresholding parameter tau by controling FDP
  if (sum(FDPprop <= alpha) > 0) tauprop = min(c(2, tau[FDPprop <= alpha]))
  if (sum(FDPprop <= alpha) == 0) tauprop = 2
  Threshold = tauprop * sqrt(kappa * log(p) / n) * (1 - EstThresh^2)
  SRec = 1 * (abs(Est) > Threshold);
  return(SRec)
}
