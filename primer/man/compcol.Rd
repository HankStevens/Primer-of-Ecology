\name{compcol}
\alias{compcol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Two-species Competition-colonization Metapopulation Model }
\description{
 This model implements a Levins-type metapopulation model for two
 species, after Hastings (1980).  For use with \code{ode} in the \code{deSolve} package.
}
\usage{
compcol(t, y, params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{ Argument for time }
  \item{y}{ A vector for population 1 and 2 }
  \item{params}{ Vector or list of parameters }
}
\value{
   Returns a list of length one, for use with \code{ode} in the \code{deSolve} package.
  \item{Component 1 }{vector of the state variables, y.}
}
\references{
Hastings, A. (1980) Disturbance, coexistence, history, and competition for
space. \emph{Theoretical Population Biology}, \bold{18}, 363--373.

 Stevens. M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R! Series. Springer.
 }
\author{  Hank Stevens <HStevens@muohio.edu>  }
\seealso{ \code{\link{levins}}, \code{\link{compcolM}}, \code{\link{succniche}}  }
\examples{
library(deSolve)
pars <- c(c1 = .3, c2 = 1, m1 = .1, m2 = .1)
pops <- c(.1,.1)
out <- ode(y=pops, t=0:20, fun=compcol, parms = pars)
matplot(out[,1], out[,-1], type='l')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }