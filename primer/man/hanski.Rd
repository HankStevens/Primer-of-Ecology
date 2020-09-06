\name{hanski}
\alias{hanski}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Core-satellite Metapopulation Model }
\description{
  A function for the core-satellite metaapopulation dynamics, for use with \code{ode} in the \code{deSolve} package.
}
\usage{
hanski(t, y, parms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{ Argument for time }
  \item{y}{ A scalar for the population variable }
  \item{parms}{ Vector or list of parameters }
}

\value{
   Returns a list of length one, for use with \code{ode} in the \code{deSolve} package.
  \item{Component 1 }{vector of the state variable (a scalar for the proportion of sites occupied).}
}
\references{ 
Hanski, I. (1982) Dynamics of regional distribution: the core and
   satellite species hypothesis. \emph{Oikos}, \bold{38}, 210--221.

Stevens, M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R!  Series. Springer.
}
\author{  Hank Stevens <HStevens@muohio.edu>  }

\seealso{ \code{\link{gotelli}}, \code{\link{hanski}},\code{\link{lande}}, \code{\link{MetaSim}}, \code{\link{clogistic}}  }

\examples{
prms <- c(ci<- 0.15, e=0.05)
out <- ode(y=.2, times=1:100, func=hanski, parms=prms )
matplot(out[,1], out[,2], type='l', ylab="p", xlab="time")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
