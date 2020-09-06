\name{moths}
\alias{moths}
\docType{data}
\title{
Moth Species Richness
}
\description{
Data set of moth species richness from Ohio and Indiana, USA.
}
\usage{data(moths)}
\format{
  A data frame with 21 observations on the following 6 variables.
  \describe{
    \item{\code{region}}{a factor with levels \code{NCT} \code{WAP} for the two collection regions.}
    \item{\code{site}}{a factor with levels associated with county names.}
    \item{\code{area}}{a numeric vector for the area of each forest fragment}
    \item{\code{spp}}{a numeric vector, the number of moth species.}
    \item{\code{lat}}{a numeric vector; latitude}
    \item{\code{long}}{a numeric vector; longitude}
  }
}

\source{
Data proided graciously by Keith Summerville and Thomas Crist; see "A Primer of Ecology with R" for specific references.
}
\references{
 K. S. Summerville and T. O. Crist. 2003. Determinants of lepidopteran community composition and species diversity in eastern deciduous forests: roles of season, ecoregion and patch size. \emph{Oikos}, 100:134-148.
 
K. S. Summerville and T. O. Crist. 2004. Contrasting effects of habitat quantity and quality on moth communities in fragmented landscapes. \emph{Ecography}, 27:3-12.
}
\examples{
data(moths)
plot(spp ~ area, data=moths)
}
\keyword{datasets}
