\name{inference}
\alias{inference}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
inference(list, alpha = 0.05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{list}{
%%     ~~Describe \code{list} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
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
function (list, alpha = 0.05) 
{
    Est = list$Est
    tscore = list$tscore
    kappa = list$kappa
    EstThresh = list$EstThresh
    n = list$n
    p = list$p
    t0 = 2
    tau = seq(0, 3.5, 0.01)
    smax = n/2
    lentau = length(tau)
    resprop = list()
    rejectprop = c()
    for (i in 1:lentau) {
        Threshold = tau[i] * sqrt(kappa * log(p)/n) * (1 - EstThresh^2)
        SRec = 1 * (abs(Est) > Threshold)
        NoNSRec = 1 * (SRec == 0)
        resprop[[i]] = which(SRec == 1, arr.ind = TRUE)
        rejectprop = c(rejectprop, max(1, (sum(SRec) - p)))
    }
    FDPprop = 2 * (p * (p - 1)) * (1 - pnorm(tau * sqrt(log(p))))/rejectprop
    FDPresprop = c()
    if (sum(FDPprop <= alpha) > 0) 
        tauprop = min(c(2, tau[FDPprop <= alpha]))
    if (sum(FDPprop <= alpha) == 0) 
        tauprop = 2
    Threshold = tauprop * sqrt(kappa * log(p)/n) * (1 - EstThresh^2)
    SRec = 1 * (abs(Est) > Threshold)
    return(SRec)
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
