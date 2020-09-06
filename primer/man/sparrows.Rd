\name{sparrows}
\alias{sparrows}
\docType{data}
\title{Song Sparrow Data Set}
\description{
  Song Sparrow (\emph{Melospiza melodia}) counts in Darrtown, OH, USA.  From  Sauer, J. R., J.E. Hines, and  J. Fallon. 2005. The North American Breeding Bird Survey, Results and Analysis 1966--2004. Version 2005.2. USGS Patuxent Wildlife Research Center, Laurel, MD.
}
\usage{data(sparrows)}
\format{
  A data frame with 36 observations on the following 3 variables.
  \describe{
    \item{\code{Year}}{a numeric vector}
    \item{\code{Count}}{a numeric vector}
    \item{\code{ObserverNumber}}{a numeric vector}
  }
}

\source{
http://www.pwrc.usgs.gov/BBS/
}

\references{ Stevens, M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R!  Series. Springer. }

\seealso{ \code{\link{PopSim}} }
\examples{
data(sparrows)
## maybe 
str(sparrows)
plot(Count ~ Year, sparrows) 
}
\keyword{datasets}
