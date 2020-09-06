\name{succniche}
\alias{succniche}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A Four-state model of Successional Dynamics }
\description{
  This is the five-state, two-species model of the succession-niche
  model, after Pacala and Rees (1998).  For use with \code{ode} in the \code{deSolve} package.
}
\usage{
succniche(t, y, params)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{t}{ Argument for the time point at integration }
  \item{y}{  A vector of length four, for states, E, M, S, and R.  }
  \item{params}{ Vector or list of parameters }
}
\value{
   Returns a list of length one, for use with \code{ode} in the \code{deSolve} package.
  \item{Component 1 }{vector of the state variable (a scalar for the proportion of sites occupied).}
}
\references{
Pacala, S.W. and Rees, M. (1998) Models suggesting field experiments to
test two hypotheses explaining successional diversity. \emph{The
  American Naturalist}, \bold{152}, 729--737.

Stevens, M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R! Series. Springer.
 }

\author{  Hank Stevens <HStevens@muohio.edu>  }
\seealso{ \code{\link{levins}}, \code{\link{compcol}}, \code{\link{compcolM}}  }
\examples{
params.suc <- c(a=7, c=0.2, g=.1, m=0.04, D=0)
t=seq(0,50,.1)
init.suc <- c(S=0, E=0.5, M=0.5, R=0.00) # Free space is implicit, F=1-(S+E+M+R).
ccg.out <- data.frame(ode(init.suc, t, succniche, params.suc))
matplot(t, ccg.out[,-1], type="l", ylab="Relative Frequency", 
         xlab="Time", ylim=c(0,1) )
legend("right", colnames(ccg.out)[5:2], lty=4:1,  bty="n")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }