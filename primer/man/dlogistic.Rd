\name{dlogistic}
\alias{dlogistic}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Discrete Logistic Growth }
\description{
 Lotka-Volterra single species discrete population growth.
}
\usage{
dlogistic(alpha = 0.01, rd = 1, N0 = 2, t = 15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{ per capita negative density dependence (a positive value will result in the usual negative effect) }
  \item{rd}{ discrete growth increment }
  \item{N0}{ initial population size }
  \item{t}{ end time point }
}
\details{
  Of the form, \deqn{N_{t+1} = N_{t} + r_d N_{t}\left(1 - \alpha N_{t}\right)}{%
  N[t+1] = N[t] + r[d] N[t] (1 - alpha N[t])}
}
\value{
 Returns a vector of population sizes from \eqn{N_0,\, \ldots,\, N_t}{N[0] to N[t]}, for integer time points.
}
\references{  
 Stevens. M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R! Series. Springer.}
\author{ Hank Stevens (HStevens@muohio.edu) }
\seealso{ \code{\link{clogistic}}, \code{\link{dlvcomp2}}, \code{\link{lvcompg}} }
\examples{
time <- 20
Ns <- dlogistic(t=20)
plot(0:time, Ns)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods}