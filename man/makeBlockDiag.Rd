\name{makeBlockDiag}
\alias{makeBlockDiag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeBlockDiag(blocksize = 4, p = 20, min.beta = 0.3, max.beta = 0.9)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{blocksize}{
%%     ~~Describe \code{blocksize} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{min.beta}{
%%     ~~Describe \code{min.beta} here~~
}
  \item{max.beta}{
%%     ~~Describe \code{max.beta} here~~
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
function (blocksize = 4, p = 20, min.beta = 0.3, max.beta = 0.9) 
{
    require(Matrix)
    require(corpcor)
    reps = p/blocksize
    S = list()
    for (i in 1:reps) {
        bd = matrix(runif(1, min.beta, max.beta), blocksize, 
            blocksize)
        diag(bd) = runif(1, 1, 1.25)
        while (!is.positive.definite(bd)) {
            diag(bd) = diag(bd) + 0.01
        }
        S[[i]] = bd
    }
    as.matrix(bdiag(S))
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
