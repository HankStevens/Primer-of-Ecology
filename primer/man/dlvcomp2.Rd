\name{dlvcomp2}
\alias{dlvcomp2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Discrete Logistic Competition }
\description{
 A function to simulate discrete 2 species Lotka-Volterra competition.
}
\usage{
dlvcomp2(N, alpha, rd = c(1, 1))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{ a vector of length two, containing N[t] for both species. }
  \item{alpha}{ a 2 x 2 matrix of the magnitudes of the per capita (negative) effects (i.e., positive value results in a negative effect). }
  \item{rd}{ a vector of length 2 containing the discrete growth increments for 2 species. }
}

\value{
  Returns a vector of length 2, containing N[t+1] for both species.
}

\references{  
 Stevens. M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R! Series. Springer.}

\author{ Hank Stevens (HStevens@muohio.edu) }

\seealso{\code{\link{dlogistic}}, \code{\link{lvcomp2}}, 
\code{\link{lvcomp3}}, \code{\link{lvcompg}} 
}

\examples{
alphs <- matrix(c( .01, .005,
                  .008, .01), ncol=2, byrow=TRUE)
t <- 20
N <- matrix(NA, nrow=t+1, ncol=2)
N[1,] <- c(10,10)
for(i in 1:t) N[i+1,] <- dlvcomp2(N[i,], alphs)
matplot(0:t, N, type='l', col=1, ylim=c(0,110))
abline(h=1/alphs[1,1], lty=3)
text(0, 1/alphs[1,1], "K", adj=c(0,0))
legend("right", c(expression("Sp.1 "*(alpha[21]==0.008)), 
                  expression("Sp.2 "*(alpha[12]==0.005))), lty=1:2, bty='n')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }
