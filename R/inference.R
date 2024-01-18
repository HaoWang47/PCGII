#' Inference() is the function to conduct simultaneous inference of estimated partial correlations.
#'
#' @param list A list returned by either `PCGII()` or `clevel()`.
#' @param alpha A pre-determined False Discovery Rate. Nominal FDR is set at 0.05 by default.
#' @returns An adjacency matrix of significant partial correlations.
#' @examples
#' PCGII_out=PCGII(df=X, prior=as.data.frame(prior_set), lambda = lam)
#' sigs_mat=inference(PCGII_out, alpha=.1)
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
