\name{compcolM}
\alias{compcolM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multi-species Competition-colonization Model, With Habitat Destruction }
\description{
  Multi-species competition colonization model, with habitat
  destruction, after Nee and May (1992).  For use with \code{ode} in the \code{deSolve} package.
}
\usage{
compcolM(t, y, params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{ Argument for each time point }
  \item{y}{  A vector for the populations }
  \item{params}{ Vector or list of parameters }
}
\value{
   Returns a list of length one, for use with \code{ode} in the \code{deSolve} package.
  \item{Component 1 }{vector of the state variables, y.}
}
\references{
Nee, S. and May, R.M. (1992) Dynamics of metapopulations: habitat destruction
and competitive coexistence. \emph{Journal of Animal Ecology},
\bold{61}, 37--40.

 Stevens. M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R! Series. Springer.
 }
\author{  Hank Stevens <HStevens@muohio.edu>  }
\seealso{ \code{\link{levins}}, \code{\link{compcol}}, \code{\link{succniche}}  }
\examples{
library(deSolve)
S <- 10
ci <- 2^seq(-5, 5, length=S)
m <- rep(.1, S)
params <- list(ci=ci, m=m, S=S, D=0)
init.N <- rep(0.01, S); t=seq(1, 200, .1)
cc.out <- ode(init.N, t, compcolM, params)
matplot(t, cc.out[, -1], type="l", ylab="Proportion of Habitat", xlab="Years")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
