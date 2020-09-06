\name{lande}
\alias{lande}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{ A Metapopulation Model with Habitat Destruction }
\description{
   A function for the metaapopulation dynamics, for use with \code{ode} in the \code{deSolve} package.
}
\usage{
lande(t, y, parms)
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
P. Kareiva and Wennergren, U. (1995) Connecting landscape patterns to ecosystem and population processes. \emph{Nature}, \bold{373}, 299--302.

Lande, R. (1987) Extinction thresholds in demographic models of territorial populations. \emph{The American Naturalist}, \bold{130}, 624--635.

Stevens, M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R!  Series. Springer.
}
\author{  Hank Stevens <HStevens@muohio.edu>  }
\seealso{ \code{\link{gotelli}}, \code{\link{hanski}},\code{\link{levins}}, \code{\link{MetaSim}}, \code{\link{clogistic}}  }
\examples{
## The function is currently defined as
function (t, y, parms) 
{
    p <- y[1]
    with(as.list(parms), {
        dp <- ci * p * (1 - D - p) - e * p
        return(list(dp))
    })
  }
library(deSolve)
p <- c(ci=.1, e=.01, D=.5)
time <- 1:10
initialN <- .3
out <- ode(y=initialN, times=time, func=lande, parms=p) 
plot(time, out[,-1], type='l')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }