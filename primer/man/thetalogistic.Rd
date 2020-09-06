\name{thetalogistic}
\alias{thetalogistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Continuous Theta-Logistic Growth }
\description{
  A function for continuous theta-logistic growth, for use with \code{ode} in the \code{deSolve} package.
}
\usage{
thetalogistic(times, y, parms)
}

\arguments{
  \item{times}{ Times points that will return N }
  \item{y}{ N }
  \item{parms}{ a vector of logistic growth parameters }
}
\details{
The user does not put these directly into this function, but rather uses \code{ode} in the \code{deSolve} package. See \code{ode} in the \code{deSolve} package.
}
\value{
  Returns of list of one component (required by \code{ode}).
}

\references{
Sibly, R.M., Barker, D., Denham, M.C., Hone, J., and Pagel, M. (2005) On
the regulation of populations of mammals, birds, fish, and
insects. \emph{Science}, \bold{309}, 607--610.

Stevens, M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R! Series. Springer.
 }
\author{ Hank Stevens (HStevens@muohio.edu) }

\seealso{\code{\link{clogistic}}}
\examples{
library(deSolve)
p <- c(r=1,alpha=.01, theta=.5)
time <- seq(1,10, by=.1)
initialN <- 10
out <- ode(y=initialN, times=time, func=thetalogistic, parms=p) 
plot(time, out[,-1], type='l')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }