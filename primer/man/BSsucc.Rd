\name{BSsucc}
\alias{BSsucc}
\docType{data}
\title{Secondary succession data}
\description{
  Functional group abundances (herbaceous annual, herbaceous perennial,
  woody) from one of the fields in the Buell-Small long term succession
  study (http://www.ecostudies.org/bss/).  Data are based on visual estimates of percent cover, using annual means of each species, which are then summed for each functional group.
}
\usage{data(BSsucc)}
\format{A data frame with 147 observations on the following 3 variables.
  \describe{
    \item{\code{AGE}}{a numeric vector; indicates the age of succession since abandonment from agriculture.}
    \item{\code{variable}}{plant functional type; a factor with levels \code{Annual}, \code{Perennial}, \code{Woody}}
    \item{\code{value}}{a numeric vector}
  }
}

\source{
 Data graciously provided by Scott Meiners (Eastern Illinois University, http://www.ecostudies.org/bss/).
}
\references{
Stevens, M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R! Series. Springer.
}
\examples{
data(BSsucc)
library(lattice)
print(xyplot(value~ AGE, groups=variable, data=BSsucc, 
             type='smooth', span=.3, ylab="Percent Cover",
             xlab="Age Since Abandonment (y)",
             auto.key=list(columns=3, space="top", lines=TRUE, points=FALSE)))
}
\keyword{datasets}
