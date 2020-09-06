\name{levins}
\alias{levins}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Classic Metapopulation Model }
\description{
  A function for the classic metaapopulation dynamics, for use with \code{ode} in the \code{deSolve} package.
}
\usage{
levins(t, y, parms)
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
Levins, R. (1969) Some demographic and genetic consequences of environmental heterogeneity for biological control. \emph{Bulletin of the Entomological Society of America}, \bold{15}, 237--240.

Stevens, M.H.H. (2009) A Primer of Ecology with R. Use R!  Series. Springer.
 }
\author{  Hank Stevens <HStevens@muohio.edu>  }
\seealso{ \code{\link{gotelli}}, \code{\link{hanski}},\code{\link{lande}}, \code{\link{MetaSim}}, \code{\link{clogistic}}  }
\examples{
## The function is currently defined as
function (t, y, parms) 
{
    p <- y[1]
    with(as.list(parms), {
        dp <- ci * p * (1 - p) - e * p
        return(list(dp))
    })
  }
library(deSolve)
p <- c(ci=.1, e=.01)
time <- 1:10
initialN <- .3
out <- ode(y=initialN, times=time, func=levins, parms=p) 
plot(time, out[,-1], type='l')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }