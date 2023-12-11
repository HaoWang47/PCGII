\name{clevel}
\alias{clevel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
clevel(df, lambda)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{df}{
%%     ~~Describe \code{df} here~~
}
  \item{lambda}{
%%     ~~Describe \code{lambda} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (df, lambda) 
{
    n = dim(df)[1]
    p = dim(df)[2]
    t0 = 2
    IndMatrix = matrix(1, p, p) - diag(rep(1, p))
    Eresidual = matrix(0, n, p)
    CoefMatrix = matrix(0, p, p - 1)
    meanX = colMeans(df)
    X = t(t(df) - meanX)
    XS = matrix(0, n, p)
    for (i in 1:p) {
        XS[, i] = X[, i]/sd(X[, i])
    }
    colnames(X) = colnames(df)
    colnames(XS) = colnames(df)
    if (missing(lambda)) {
        shat = sqrt(n/(log(p)^3))
        lambda = sqrt(2 * (2 + 0.01) * log(p/shat)/n)
    }
    for (i in 1:p) {
        out = glmnet(XS[, -i], X[, i], family = "gaussian", lambda = lambda)
        Coef = out$beta
        Predict = predict(out, XS[, -i], type = "link")
        CoefMatrix[i, ] = as.numeric(Coef)/apply(X[, -i], 2, 
            sd)
        Eresidual[, i] = X[, i] - Predict
    }
    CovRes = t(Eresidual) \%*\% Eresidual/n
    Est = matrix(1, p, p)
    for (i in 1:(p - 1)) {
        for (j in (i + 1):p) {
            temp = Eresidual[, i] * Eresidual[, j] + Eresidual[, 
                i]^2 * CoefMatrix[j, i] + Eresidual[, j]^2 * 
                CoefMatrix[i, j - 1]
            Est[i, j] = mean(temp)/sqrt(diag(CovRes)[i] * diag(CovRes)[j])
            Est[j, i] = Est[i, j]
        }
    }
    EstThresh = Est * (abs(Est) >= (t0 * sqrt(log(p)/n) * IndMatrix))
    kappa = (n/3) * mean(colSums(Eresidual^4)/(colSums(Eresidual^2))^2)
    SE = sqrt((kappa * (1 - EstThresh^2))^2/n)
    tscore = Est/SE
    return(list(Est = Est, tscore = tscore, kappa = kappa, EstThresh = EstThresh, 
        n = n, p = p))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
