\name{unMat}
\alias{unMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
unMat(X_est, X_p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X_est}{
%%     ~~Describe \code{X_est} here~~
}
  \item{X_p}{
%%     ~~Describe \code{X_p} here~~
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
function (X_est, X_p) 
{
    p = dim(X_est)[2]
    out = matrix(NA, nrow = p * (p - 1)/2, ncol = 4)
    ind = 1
    for (i in 1:(p - 1)) {
        for (j in (i + 1):p) {
            out[ind, 1] = i
            out[ind, 2] = j
            out[ind, 3] = X_est[i, j]
            out[ind, 4] = X_p[i, j]
            ind = ind + 1
        }
    }
    colnames(out) = c("node1", "node2", "pcor", "pval")
    out
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
