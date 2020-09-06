\name{MetaSim}
\alias{MetaSim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulation of  Stochastic Metapopulation Models }
\description{
  Originally focused on creating a community of core-satellite species,
  this function allows simulation of several metapopulation models,
  where colonization and extinction rates are stochastic draws from
  uniform distributions, with specified means and ranges.
}
\usage{
MetaSim(Time = 50, NSims = 1, method = "hanski",
               ci = 0.25, e = 0.25, phi = 0.75, p0 = 0.5, D = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Time}{ A scalar for the number of time steps over which to simulate each population. }
  \item{NSims}{ A scalar for the number of simulations, which is analogous to the number of species in the community. }
  \item{method}{ A character string, in quotes, specifying which metapopulation model to use: "hanksi", "gotelli", "lande", "levins". See "See Also" below. }
  \item{ci}{ Scalar for mean colonization rate. }
  \item{e}{ Scalar for mean extinction rate. }
  \item{phi}{ A scalar for the relative variability in rates. See Details. }
  \item{p0}{ Initial proportion of sites occupied for each species. }
  \item{D}{ Parameter for habitat destruction; applies to only the "lande" model. }
}
\details{
  \code{phi} is one half of the relative range of each rate
  (colonization and extinction).  For each time step, each rate is drawn
  from a uniform distribution, Unif(rate-phi, rate+phi).  Thus, the range is 2*phi, and center on the specified mean (ci or e).
}
\value{
Function returns a list with these components.
  \item{method }{The method used (default is "hanski").}
  \item{time }{The integer sequence of times, from 0 to the value of the argument \code{Time}.}
\item{Ns }{\code{Time} by \code{NSims} matrix of observed population sizes.}
\item{Parameters }{A named vector of the parameters used for the simulations.} }

\references{ 
Stevens. M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R! Series. Springer.
}
\author{  Hank Stevens <HStevens@muohio.edu>  }
\seealso{ \code{\link{gotelli}}, \code{\link{hanski}}, \code{\link{levins}} }
\examples{
out <- MetaSim(NSims=2)
pops <- out$Ns
matplot(out$t, pops, type='l')
title(sub=paste(out$method, "model"))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }