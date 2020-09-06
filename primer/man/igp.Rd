\name{igp}
\alias{igp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A Lotka-Volterra Model of Intraguild Predation }
\description{
 A Lotka-Volterra model of intraguild predation, after Holt and Polis
 (1997).  For use with \code{ode} in the \code{deSolve} package.
}
\usage{
igp(t, y, params)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{t}{ the time point for a realization of the integration.}
  \item{y}{ the vector of populations, at each time t. }
  \item{params}{ a vector or list containing the necessary parameters. }
}

\value{
 Returns a list of length one which is the vector of the rates of increase (required by \code{ode}).
}
\references{ 
Holt, R.D. and  Polis, G.A. (1997) A theoretical framework for
 intraguild predation. \emph{The American Naturalist}, \bold{149}, 745--764.

Stevens, M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R! Series. Springer.
 }
\author{ Hank Stevens <HStevens@muohio.edu> }
\seealso{\code{\link{lvcompg}}, \code{\link{scheffer}} }
\examples{
library(deSolve)
params <- c(bpb= 0.032, abp=10^-8, bpn=10^-5, anp=10^-4,  mp=1,
             bnb=0.04, abn=10^-8, mn=1,
             r=1, abb=10^-9.5)
t=seq(0, 60, by=.1)
N.init <- c(B = 10^9, N = 10^4, P = 10^3)
igp.out <- ode(N.init, t, igp, params)
matplot(t, log10(igp.out[,-1]+1), type="l", 
          ylab="log(Abundance)")
legend('right', c("B", "N", "P"), lty=1:3, col=1:3, lwd=2,
       bty="n")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }