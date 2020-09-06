\name{ClostExp}
\alias{ClostExp}
\docType{data}
\title{Closterium Population Data}
\description{
   Data set from an experiment testing whether trophic heterogeneity influences long term population sizes. 
}
\usage{data(ClostExp)}
\format{
  A data frame with 144 observations on the following 5 variables.
  \describe{
    \item{\code{Nutrients}}{a factor with levels \code{high} \code{low}}
    \item{\code{No.per.ml}}{a numeric vector}
    \item{\code{Day}}{a numeric vector}
    \item{\code{rep}}{a factor with levels \code{a} \code{b} \code{c} \code{d}}
    \item{\code{ID}}{an ordered factor with levels \code{a.low} < \code{d.low} < \code{c.low} < \code{b.low} < \code{c.high} < \code{a.high} < \code{d.high} < \code{b.high}}
  }
}
\source{
  Stevens, M. H. H. and Steiner, C. E. (2006) Effects of predation and
  nutrient enrichment on a food web with edible and inedible
  prey. \emph{Freshwater Biology}, \bold{51}, 666--671.
}
\references{
   Stevens, M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R! Series. Springer.
}
\examples{
data(ClostExp)
library(lattice)
xyplot(log10(No.per.ml) ~ Day|Nutrients, data=ClostExp, groups=rep, type='l')
}
\keyword{datasets}
