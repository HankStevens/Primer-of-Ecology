\name{lvcomp2}
\alias{lvcomp2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Two Species Lotka-Volterra Competition }
\description{
 System of ordinary differential equations for two species Lotka-Volterra competition. For use with \code{ode} in the \code{deSolve} package.
}
\usage{
lvcomp2(t, n, parms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{ Times points that will return N. }
  \item{n}{a vector of length two for the population sizes at time  = t. }
  \item{parms}{ vector or list of model parameters (see details below). }
}
\details{
 The parameters include \code{r1, r2, a11, a12, a21, a22} with the usual meanings. Here the a's are the per capita effects which determine K (\code{a11 = 1/K}).
}
\value{
 Returns a list of length one which is the rate of increase (required by \code{ode}).
}
\references{
Lotka, A.J. (1956) \emph{Elements of Mathematical Biology}. Dover Publications, Inc.

Stevens. M.H.H. (2009) A Primer of Theoretical Population Ecology with R. Use R! Series. Springer.
 }

\author{  Hank Stevens <HStevens@muohio.edu>  }

\seealso{ \code{\link{lvcomp3}}, \code{\link{lvcompg}}, \code{\link{clogistic}} }
\examples{
## The function is currently defined as
function (t, n, parms) 
{
    with(as.list(parms), {
        dn1dt <- r1 * n[1] * (1 - a11 * n[1] - a12 * n[2])
        dn2dt <- r2 * n[2] * (1 - a22 * n[2] - a21 * n[1])
        list(c(dn1dt, dn2dt))
    })
  }
library(deSolve) 
parms <- c(r1 = 1, r2 = 0.1, a11 = 0.2, a21 = 0.1, a22 = 0.02, a12 = 0.01) 
initialN <- c(1, 1) 
out <- ode(y = initialN, times = 1:100, func = lvcomp2, parms = parms) 
matplot(out[, 1], out[, -1], type = "l") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }

