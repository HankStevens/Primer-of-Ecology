\name{lvcomp3}
\alias{lvcomp3}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Three Species Lotka-Volterra Competition }
\description{
 System of ordinary differential equations for three species Lotka-Volterra competition. For use with \code{ode} in the \code{deSolve} package.
}
\usage{
lvcomp3(t, n, parms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{ the time for each integration. }
  \item{n}{ a vector of length three for the population sizes at time  = t. }
  \item{parms}{ vector or list of model parameters (see details below). }
}

\details{
 The parameters include \code{r, a} with the usual meanings. Here the \code{a}'s are the per capita effects which determine K (\code{a11 = 1/K1}).
}
\value{
 Returns a list of length one which is the rate of increase (required by \code{ode}).
}

\references{
Lotka, A.J. (1956) \emph{Elements of Mathematical Biology}. Dover Publications, Inc.

Stevens. M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R! Series. Springer.
 }

\author{  Hank Stevens <HStevens@muohio.edu>  }

\seealso{ \code{\link{lvcomp2}}, \code{\link{lvcompg}}, \code{\link{clogistic}} }
\examples{
## The function is currently defined as
function (t, n, parms) 
{
    with(as.list(parms), {
        dn1dt <- r1 * n[1] * (1 - a11 * n[1] - a12 * n[2] - a13 * 
            n[3])
        dn2dt <- r2 * n[2] * (1 - a22 * n[2] - a21 * n[1] - a23 * 
            n[3])
        dn3dt <- r3 * n[3] * (1 - a33 * n[3] - a31 * n[1] - a32 * 
            n[2])
        list(c(dn1dt, dn2dt, dn3dt))
    })
  }

library(deSolve) 
parms <- c(r1 = 0.1, r2 = 0.2, r3 = 0.3, 
a11 = 0.1, a12 = 0.01, a13 = 0.01, 
a21 = 0.01, a22 = 0.15, a23 = 0.01, 
a31 = 0.01, a32 = 0.01, a33 = 0.2) 
initialN <- c(1, 1, 1) 
out <- ode(y = initialN, times = 1:100, func = lvcomp3, parms = parms) 
matplot(out[, 1], out[, -1], type = "l") 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
