\name{LVCompGames}
\alias{LVCompGames}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Lotka-Volterra Competition Games }
\description{
 A pedagogical tool for learning about isoclines and stable equilibria. 
}
\usage{
LVCompGames(Alpha=NULL, r1=.1, r2=.1, num=20, time=10, step=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Alpha}{a 2 x 2 matrix of coefficients; if NULL, then a sensible
    random matrix will be generated - fun for games.}
  \item{r1}{ the intrinsic rate of increase for species 1 (isocline -- a
    solid line). }
  \item{r2}{ the intrinsic rate of increase for species 2 (isocline -- a
    dashed line). }
  \item{num}{ the desired number of random starting points for trajectories. }
  \item{time}{ the number of time steps for each trajectory. }
  \item{step}{ the number of time steps for which the integration is
    estimated (has no effect on accuracy, merely the smoothness of the line). }
}

\value{
First generates a plot of isoclines and initial points; the user is then
prompted to "Hit <return>" in the console to see the trajectories. These
will indicate the rate and path the trajectories, demonstrating, among
other things, whether the equilibrium is stable.
}
\references{
Lotka,  A.J. (1956) \emph{Elements of Mathematical Biology}. Dover Publications, Inc.

Stevens, M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R! Series. Springer.
 }

\author{  Hank Stevens <HStevens@muohio.edu>  }
\seealso{\code{\link{lvcomp2}}, \code{\link{lvcompg}}, \code{\link{clogistic}}, }
\examples{
LVCompGames() # Hit return in the console to see the trajectories.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }