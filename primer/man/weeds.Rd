\name{weeds}
\alias{weeds}
\docType{data}
\title{Percent cover of six perennial herbaceous plants}
\description{
 Percent cover of six of the most common herbaceous perennial species
 from the Buell-Small long term succession study
 (http://www.ecostudies.org/bss/).  Data are plot-level visual estimates of percent cover.
}
\usage{data(weeds)}
\format{
  A data frame with 15140 observations on the following 8 variables.
  \describe{
    \item{\code{X}}{a numeric vector}
    \item{\code{Genus}}{a factor with levels \code{Aster}, \code{Euthamia}, \code{Solidago}}
    \item{\code{Epithet}}{a factor with levels \code{canadensis}, \code{gigantea}, \code{graminifolia}, \code{novae-angliae}, \code{pilosus}, \code{rugosa}}
    \item{\code{FieldName}}{a factor with levels \code{C3,} \code{C4}, \code{C5}, \code{C6}, \code{C7}, \code{D1}, \code{D2}, \code{D3}, \code{E1}, \code{E2}}
    \item{\code{Age}}{a numeric vector indicating the number of years of succession, since abandonment from agriculture.}
    \item{\code{PlotId}}{a numeric vector}
    \item{\code{Cover}}{a numeric vector; percent cover, estimated visually.}
    \item{\code{Species}}{a factor with levels \code{A.novae-angliae}, \code{A.pilosus}, \code{E.graminifolia}, \code{S.canadensis}, \code{S.gigantea}, \code{S.rugosa}}
  }
}

\source{
Data graciously provided by Scott Meiners (Eastern Illinois University,
http://www.ecostudies.org/bss/).
}
\references{
Stevens, M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R! Series. Springer.
}
\examples{
data(weeds)
str(weeds); 
library(lattice)
xyplot(Cover ~ Age, data=weeds, groups=Species, type=c("a"),
auto.key=list(lines=TRUE, points=FALSE, columns=3),
ylim=c(-1,20))
}
\keyword{datasets}
