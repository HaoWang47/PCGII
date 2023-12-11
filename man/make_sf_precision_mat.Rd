\name{make_sf_precision_mat}
\alias{make_sf_precision_mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
make_sf_precision_mat(e = 1, p = 20, lower = 0.2, upper = 0.5, diag = 0.1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{e}{
%%     ~~Describe \code{e} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{lower}{
%%     ~~Describe \code{lower} here~~
}
  \item{upper}{
%%     ~~Describe \code{upper} here~~
}
  \item{diag}{
%%     ~~Describe \code{diag} here~~
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
function (e = 1, p = 20, lower = 0.2, upper = 0.5, diag = 0.1) 
{
    require(igraph)
    g <- sample_pa(p, power = 1, m = e, directed = FALSE)
    omega = as_adjacency_matrix(g) \%>\% as.matrix()
    for (h1 in 1:(p - 1)) {
        for (h2 in (h1 + 1):p) {
            if (omega[h1, h2] != 0) {
                temp = runif(1, lower, upper) * sample(c(-1, 
                  1), size = 1)
                omega[h1, h2] = temp
                omega[h2, h1] = temp
            }
        }
    }
    diag(omega) = rowSums(abs(omega)) + diag
    omega
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
