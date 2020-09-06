\name{CandG}
\alias{CandG}
\docType{data}
\title{ Data drawn approximately from Collins and Glenn (1991)}
\description{
 Numbers of species which were observed in 1--19 sites.
}
\usage{data(CandG)}
\format{
  The format is:
 num [1:19] 32 16 10 9 8 7 8 6 4 5 ...
}
\source{Data are approximations of histogram values for a single year of data in:

Collins, S.L. and Glenn, S.M. (1991) Importance of spatial and temporal
dynamics in species regional abundance and distribution. \emph{Ecology},
\bold{72}, 654--664.

}
\references{
Stevens, M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R! Series. Springer.
}
\examples{
data(CandG)
barplot(CandG, names=1:19)
}
\keyword{datasets}
