\name{ProjMat}
\alias{ProjMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Palm Data Resampling }
\description{
  A function to create a projection matrix from the Chamaedorea palm data sets. Specific to the forms of these data sets.
}
\usage{
ProjMat(stagedat, fruitdat, seeddat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{stagedat}{ Data on the vegetative stage of individuals through time; see \code{\link{stagedat}}. }
  \item{fruitdat}{ Data on fruit production of individuals; see \code{\link{fruitdat}}. }
  \item{seeddat}{ Data on seed fates; see \code{\link{seeddat}}. }
}

\value{
 Returns a 5 x 5 projection matrix based on palm data.
}
\references{
Caswell, H. (2001) Matrix Population Models: Construction, Analysis, and Interpretation. Sinauer Associates, Inc., Sunderland, MA, USA, 2nd ed.

Ellner, S.P. and Guckenheimer, J. (2006) Dynamic Models in Biology. Princeton University Press.

Stevens, M.H.H. (2009) A Primer of Ecology with R. Use R! Series. Springer.

Berry, E.J., Gorchov, D.L., Endress, B.A.,  and Stevens, M.H.H. (2008) Source-sink dynamics within a plant population: the impact of substrate and herbivory on palm demography. Population Ecology, \bold{50}, 63--77.
 }


\seealso{\code{\link{DemoInfo}}, \code{\link{DemoBoot}},
		\code{\link{stagedat}},
		\code{\link{fruitdat}},
		\code{\link{seeddat}},
  }
\examples{
data(stagedat); data(fruitdat); data(seeddat)
ProjMat(stagedat, fruitdat, seeddat)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }