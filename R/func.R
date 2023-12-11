# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Cmd + Shift + B'
#   Check Package:             'Cmd + Shift + E'
#   Test Package:              'Cmd + Shift + T'

# PCGII
##########################################################################################
##########################################################################################
# This script contains all the R functions needed to implement the approach PCGII
# Technical details please see the manuscript.
##########################################################################################
##########################################################################################



require(tidyverse)

sigs2vec=function(sigs, P){
  # A utility function
  # This function takes PCGII inference results as input, generates a
  # adjacecy matrix corresponding to the significant partial correlations
  # and returns a vector of all off-diagonal elements. The function does
  # require package corpcor.
  require(corpcor)
  m=matrix(0,P,P)
  for (h in 1: dim(sigs)[1]){
    m[sigs[h,1],sigs[h,2]]=1
  }
  sm2vec(m)
}

sigs2mat=function(sigs, P){
  # A utility function
  # This function takes PCGII inference results as input, generates a
  # adjacecy matrix corresponding to the significant partial correlations
  # and returns the matrix. The function does require package corpcor.
  require(corpcor)
  m=matrix(0,P,P)
  for (h in 1: dim(sigs)[1]){
    m[sigs[h,1],sigs[h,2]]=1
    m[sigs[h,2],sigs[h,1]]=1
  }
  m
}

unMat=function(X_est, X_p){
  # This function takes the estimated p by p partial correlations matrix and
  # corresponding element-wise pvalue matrix as input, and return a matrix where
  # each row corresponds to a pair of nodes with the estimated partial correlation
  # and p-value.

  # X is a p by p matrix
  p=dim(X_est)[2]
  out=matrix(NA, nrow=p*(p-1)/2, ncol=4)
  ind=1
  for (i in 1:(p-1)){
    for (j in (i+1):p){
      out[ind,1]=i
      out[ind,2]=j
      out[ind,3]=X_est[i,j]
      out[ind,4]=X_p[i,j]
      ind=ind+1
    }
  }
  colnames(out)=c("node1","node2","pcor","pval")
  out
}

make_sf_precision_mat=function(e=1, p=20, lower=.2, upper=.5, diag=0.1){
  require(igraph)
  g <- sample_pa(p, power=1, m=e, directed = FALSE)
  omega=as_adjacency_matrix(g) %>% as.matrix()
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

make_random_precision_mat=function(eta=.01, p=20, lower=.2, upper=.5, diag=0.1){
  require(igraph)
  g <- sample_gnp(p, eta, directed = FALSE)
  omega=as_adjacency_matrix(g) %>% as.matrix()
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

makeBlockDiag=function(blocksize=4, p=20, min.beta=0.3, max.beta=0.9){ # blocksize has to be a factor of p
  require(Matrix)
  require(corpcor)
  reps=p/blocksize
  S=list()
  for (i in 1:reps) {
    bd=matrix(runif(1,min.beta, max.beta), blocksize, blocksize)
    diag(bd)=runif(1,1,1.25)
    while(!is.positive.definite(bd)){diag(bd)=diag(bd)+0.01}
    S[[i]]=bd
  }
  as.matrix(bdiag(S))
}

undirected_prior = function(prior){
  ## An utility function to pre-process the input prior set. This function will ensure the input prior set corresponds to an undirected prior network. If the prior network is believed to be directed, no pre-processing of the prior set is needed.

  ## Input:
  #   - prior: the prior set, a k by 2 matrix, in which each row corresponds to a pair of nodes (any omics features) that are connected under prior belief.

  ## Output: This function returns
  ##  - a pre-processed prior set, in which the connection between any pair of nodes is undirected.

  ## Remark: this function is not necessary. Prior set should be considered carefully before running the network analysis. If the prior network connections are believed to be undirected while the prior set only includes one way connections for simplicity, this function will duplicate the connections and swap the direction automactically.
  rbind.data.frame(prior, prior %>% transform(row = pmax(row, col), col = pmin(row, col))) %>%
    arrange(row, col) %>% unique()
}






## PCGII() is the function to apply the proposed method to get the estimated partial correlation graph with information incorporation

## Input:
#   - df: the main expression dataset, an n by p matrix, in which each row corresponds to a sample and each column represents expression/abundance of an omics feature.
#   - prior: the prior set, a k by 2 dataframe, in which each row corresponds to a pair of nodes (any omics features) that are connected under prior belief. Note, prior input has to be dataframe.
#   - lambda: the regularization parameter, used in the node-wise regression. If missing, default lambda will be used which is at the order of sqrt(2*log(p)/n).

## Output: This function returns a list of
#  - estimated partial correlation matrix (Est),
#  - sparse partial correlation estimation matrix with threshold (EstThresh),
#  - estimated kappa (kappa),
#  - estimated test statistics matrix of partial correlations (tscore),
#  - sample size (n) and number of nodes (p).

## Remark: mathematical standardization will be automatically done within the function.
PCGII=function(df, prior, lambda){
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

## clevel() is the function to apply the method originally proposed in paper '> Qiu, Y., & Zhou, X. H. (2020). Estimating c-level partial correlation graphs with application to brain imaging. Biostatistics (Oxford, England), 21(4), 641–658. https://doi.org/10.1093/biostatistics/kxy076', code credited to Dr. Yumou Qiu.

## Input:
#   - df: the main expression data, an n by p matrix, in which each row corresponds to a sample and each column represents expression/abundance of an omics feature.
#   - lambda: the regularization parameter, used in the node-wise regression. If missing, default lambda will be used which is at the order of sqrt(2*log(p)/n).

## Output: This function returns a list of
#   - estimated partial correlation matrix (Est),
#   - sparse partial correlation estimation matrix with threshold (EstThresh), e
#   - stimated kappa (kappa),
#   - estimated test statistics matrix of partial correlations (tscore),
#   - sample size (n) and number of nodes (p).

## Remark: mathematical standardization will be automatically done within the function.
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

## Inference() is the function to conduct simultaneous inference of estimated partial correlations. The detail of this inference approach can be found in the manuscript or in 'Qiu, Y., & Zhou, X. H. (2020). Estimating c-level partial correlation graphs with application to brain imaging. Biostatistics (Oxford, England), 21(4), 641–658. https://doi.org/10.1093/biostatistics/kxy076'. Code credited to Dr. Yumou Qiu

## Input:
#   - list: a list returned by either `PCGII()` or `clevel()`.
#   - alpha: pre-determined False Discovery Rate. Nominal FDR is set at 0.05 by default.

## Output:
# a list contains the dataframe of pairs with significant partial correlations.
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


# unique_cord=function(cord){
#   cord %>% as.data.frame %>% transform(row = pmin(row, col), col = pmax(row, col)) %>%
#     arrange(row, col) %>%
#     unique()
# }
