\name{PopSim}
\alias{PopSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Population Simulator }
\description{
  Resampling stochastic simulator for a single density-independent population.
}
\usage{
PopSim(Rs, N0, years = 50, sims = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Rs}{ vector of observed annual growth rates (N[t+1]/N[t]). }
  \item{N0}{ initial population size. }
  \item{years}{ number of years to simulate. }
  \item{sims}{ number of simulated populations. }
}
\details{
  Designed to simulate trajectories based on resampled observed N[t+1]/N[t].
}
\value{
 Returns a matrix of population sizes for time = t (rows) for each replicated population (columns).
}
\references{Stevens,  M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R!  Series. Springer. }

\seealso{ \code{\link{sparrows}} }
\examples{
lambdas <- runif(10, .5, 1.5)
out <- PopSim(Rs=lambdas, years=50, N0=50)
matplot(0:50, out, type='l')
summary(out[51,])
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }