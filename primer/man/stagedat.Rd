\name{stagedat}
\alias{stagedat}
\docType{data}
\title{Chamaedorea Palm Vegetative Stage Data Set}
\description{
 These stage data provide the stage of each individual in the study. 
Each row is 
an individual, and its ID number is in column 1. Data in columns 2-4 identify 
its stage in years 2003-2005. 
}
\usage{data(stagedat)}
\format{
  A data frame with 414 observations (414 individual plants) on the following 4 variables.
  \describe{
    \item{\code{PalmNo}}{a numeric vector; plant ID number}
    \item{\code{Y2003}}{a numeric vector; stage of individual i in 2003}
    \item{\code{Y2004}}{a numeric vector; stage of individual i in 2004}
    \item{\code{Y2005}}{a numeric vector; stage of individual i in 2005}
  }
}

\source{
Berry, E.J.,Gorchov, D.L., Endress, B.A., and Stevens, M.H.H. (2008) Source-sink dynamics within a plant population: the impact of substrate and herbivory on palm demography. \emph{Population Ecology}, \bold{50}, 63--77.
}
\references{  
 Stevens. M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R! Series. Springer.}
\author{ Hank Stevens (HStevens@muohio.edu) }
\seealso{ \code{\link{DemoBoot}}, \code{\link{DemoInfo}}, \code{\link{seeddat}}, \code{\link{fruitdat}} }
\examples{
data(stagedat)
table(stagedat$Y2004)
}
\keyword{datasets}
