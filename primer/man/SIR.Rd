\name{SIR}
\alias{SIR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ The S-I-R Epidemilogical Disease Model }
\description{
  The S-I-R epidemilogical disease model, for use with  \code{ode} in the \code{deSolve} package.
}
\usage{
SIR(t, y, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{ times points for which values will be returned }
  \item{y}{ the vector of disease states of hosts (S, I, R) }
  \item{p}{ a vector of parameters }
}
\details{
The user does not put these directly into this function, but rather uses \code{ode} in the \code{deSolve} package.
}
\value{
  Returns of list of one component (required by \code{ode}).
}
\references{
Ellner, S.P. and Guckenheimer, J. (2006) \emph{Dynamic Models in Biology}, Princeton University Press.

Kermack, W.O. and McCormick, W.G. (1927) A contribution to the
mathematical theory of epidemics. \emph{Proceedings of the Royal
  Society, Series A}, \bold{115}, 700--721.

Stevens, M.H.H. (2009) \emph{A Primer of Ecology with R}, Use R! Series. Springer.
 }
\author{  Hank Stevens <HStevens@muohio.edu>  }
\seealso{ \code{\link{ross}}, \code{\link{SIRf}}, \code{\link{SIRbd}}  }
\examples{
N <- 10^3; I <- R <- 1; S <- N - I - R
parms <- c(B=.01, g=4)
months <- seq(0, 3, by=0.01)
require(deSolve)
SIR.out <- data.frame( ode(c(S,I,R), months, SIR, parms) )
matplot(months, SIR.out[,-1], type='l')
legend('right', c('R', 'I', 'S'), lty=3:1, col=3:1, bty='n')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
