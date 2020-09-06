\name{clogistic}
\alias{clogistic}

\title{ Continuous Logistic Growth }
\description{
  A function for continuous logistic growth, for use with \code{ode} in the \code{deSolve} package.
}
\usage{
clogistic(times, y, parms)
}

\arguments{
  \item{times}{ Times points that will return N }
  \item{y}{ N }
  \item{parms}{ a vector of logistic growth parameters }
}
\details{
The user does not put these directly into this function, but rather uses
\code{ode} in the \code{deSolve} package.

The function is based on the logistic growth equation
\deqn{\frac{dN}{dt} = rN\left(1-\alpha N\right)}{%
  dN/dt = rN(1-alpha N)}

or equivalently
\deqn{\frac{dN}{dt} = rN\left(\frac{K-N}{K}\right)}{%
  dN/dt = rN(K-N)/K}
}
\value{
  Returns of list of one component (required by \code{ode}).
}
\references{
Lotka, A.J. (1956) \emph{Elements of Mathematical Biology}. Dover Publications, Inc.

 Stevens. M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R! Series. Springer.
 }

\author{ Hank Stevens <HStevens@muohio.edu> }

\seealso{\code{\link{dlogistic}}, \code{\link{thetalogistic}}, \code{\link{levins}}, \code{\link{lvcomp2}}, 
\code{\link{lvcompg}}
}
\examples{
library(deSolve)
p <- c(1,.01)
time <- 1:10
initialN <- 10
out <- ode(y=initialN, times=time,
             func=clogistic, parms=p) 
plot(time, out[,-1], type='l')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }