\name{predpreyLV}
\alias{predpreyLV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The Lotka-Volterra Predator-prey Model }
\description{
  The Lotka-Volterra predator-prey model, for use with \code{ode} in the \code{deSolve} package.
}
\usage{
predpreyLV(t, y, params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{ Argument for time }
  \item{y}{ A vector of length 2, for population 1 and 2 }
  \item{params}{ Vector or list of parameters }
}
\value{
   Returns a list of length one, for use with \code{ode} in the \code{deSolve} package.
  \item{Component 1 }{vector of the state variables, y.}
}
\references{
Lotka, A.J. (1956) \emph{Elements of Mathematical Biology}. Dover Publications, Inc.

Stevens. M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R! Series. Springer.
 }
\author{  Hank Stevens <HStevens@muohio.edu>  }
\seealso{ \code{\link{predpreyRM}}  }
\examples{
params1 <- c(b=.5, a=.01, s=.2, e=.1) 
Time <- seq(0,100, by=.1) # Set time here
LV.out <- ode(c(H0=25,P0=5), Time, predpreyLV, params1)
matplot(Time, (LV.out[,2:3]), type="l", ylab="Population Size") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }