\name{seeddat}
\alias{seeddat}
\docType{data}
\title{ Chamaedorea Palm Seed Data Set}
\description{
This data set provides the number of seeds per individual in 2004.
These data are the fates of each seed in a sample of 400 seeds.  Seeds may have germinated (2), remained viable (1), or died (0). 
}
\usage{data(seeddat)}
\format{
  A data frame with 400 observations on the following variable.
  \describe{
    \item{\code{Y2004}}{a numeric vector}
  }
}

\source{
Berry, E.J., Gorchov, D.L., Endress, B.A.,  and Stevens, M.H.H. (2008)
Source-sink dynamics within a plant population: the impact of substrate
and herbivory on palm demography. \emph{Population Ecology}, \bold{50}, 63--77.
}
\references{  
 Stevens. M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R! Series. Springer.}
\author{ Hank Stevens (HStevens@muohio.edu) }
\seealso{ \code{\link{DemoBoot}}, \code{\link{DemoInfo}}, \code{\link{seeddat}}, \code{\link{stagedat}} }

\examples{
data(seeddat)
table(seeddat)
}
\keyword{datasets}
