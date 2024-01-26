# PCGII [![](https://img.shields.io/badge/Release-v1.1.2-blue.svg)](https://github.com/haowang47/PCGII/commits/main) [![License: MIT v3](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/license/mit/)

R Package for Information-incorporated Gene Network Construction with FDR Control

### Authors:
> Hao Wang, Yumou Qiu and Peng Liu.

### Contact:
> [haydo.wang@outlook.com] (Hao Wang)

### Citation:
> Wang, H., Qiu, Y.\*, Guo, H., Yin, Y., Liu, P.\*, 2024. Information-incorporated Gene Network Construction with FDR Control. Under review.

# Installation and Package loading
```r
# R version is required >= 3.4.4
# When the first time to use the package, please make sure dependent packages are installed under your R environment, if not, please use commands below to install
> #install.packages("tidyverse")
> #install.packages("glmnet")
> #install.packages("mvtnorm")
> #install.packages("igraph")
> #install.packages("Matrix")
# install "devtools" package in your R environment
> # devtools::install_github("HaoWang47/PCGII")
> library(PCGII)
> library(corpcor)
> library(glmnet)
> library(igraph)
> library(Matrix)
> library(mvtnorm)
> library(tidyverse)
```

This is a tutorial script for researchers who are interested in applying PCGII on omics data to learn the direct association structure of omics features. The main function `PGCII()` takes a biologically pre-processed expression data matrix as input, and returns a list of statistics (estimates and test statistics). The function `inference()` takes a list returned by `PGCII()` as input and conduct simultaneous test to identify significant partial correlations with False Discovery Rate (FDR) controlled at a pre-determined nominal level (0.05 by default).

### Usage

`PCGII()`

  - Input:
    - `df`: the main expression data, an $n$ by $p$ matrix/dataframe, in which each row corresponds to a sample and each column represents expression/abundance of an omics feature;
    - `prior`: the prior set, a $k$ by $2$ dataframe, in which each row corresponds to a pair of nodes (any omics features) that are connected under prior belief. Note, prior input has to be dataframe with column names **"row"** and **"col"**;
    - `lambda`: the regularization parameter, used in the node-wise regression. If missing, default lambda will be used which is at the order of sqrt(2log(p/sqrt(n))/n).
  - Remark: mathematical standardization will be automatically done within the function.
  - Output: This function returns a list of
    - `Est`: estimated partial correlation matrix;
    - `EstThresh`: sparse partial correlation estimation matrix with threshold;
    - `kappa`: estimated ratio of forth and squared second moment of residuals, please refer to the manuscript for details;
    - `tscore`: estimated test statistics matrix of partial correlations;
    - `n`: sample size;
    - `p`: number of genes under study.

`Inference()`

  - Input:
    - `list`: a list returned by either `PCGII()` or `clevel()`.
    - `alpha`: pre-determined False Discovery Rate. Nominal FDR is set at 0.05 by default.
  - Output: an adjacency matrix of significant partial correlations.

# Network Analysis

Simulate data $X$ from a scale-free network $g$.

```r
> # Simulating data
> set.seed(1234567)
> n=50 # sample size
> p=30 # number of nodes
>
> omega=make_random_precision_mat(eta=.01, p=p)
>
> Sigma=solve(omega) # population covariance matrix, which is used to generate data
> X = rmvnorm(n = n, sigma = Sigma) # simulate expression data
```

Network analysis of data matrix `X`.

```
> # determine tuning parameter: fixed lambda
> lam=2*sqrt(log(p)/n)
>
> # create prior set: directed prior network
> prior_set=matrix(data=c(6,5, 28,14), nrow=2, ncol=2, byrow = TRUE)
> colnames(prior_set)=c("row", "col")
> PCGII_out=PCGII(df=X, prior=as.data.frame(prior_set), lambda = lam)
> inference_out=inference(list=PCGII_out)
> diag(inference_out)=0
> # Visualization
> inference_out %>%
+   graph_from_adjacency_matrix(mode = "undirected") %>%
+   plot(vertex.size=4, vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.5)
>
> # create prior set: undirected prior network
> PCGII_out=PCGII(df=X, prior=undirected_prior(prior_set), lambda = lam)
> inference_out=inference(list=PCGII_out)
> diag(inference_out)=0
> # Visualization
> inference_out %>%
+   graph_from_adjacency_matrix(mode = "undirected") %>%
+   plot(vertex.size=4, vertex.label.dist=0.5, vertex.color="red", edge.arrow.size=0.5)
```

