\name{scheffer}
\alias{scheffer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ A resource-based model of alternative stable states }
\description{
A model of floating vs. submerged plant dominance in shallow aquatic
systems, after Scheffeer \emph{et al}. (2003).  For use with \code{ode}
in the \code{deSolve} package.  Floating plants are better competitors for light, as long as submerged plants cannot drive down nitrogen levels.
}
\usage{
scheffer(t, y, p)
}

\arguments{
  \item{t}{ the time point for a realization of the integration.}
  \item{y}{ the vector of populations, at each time t. }
  \item{p}{ a vector or list containing the necessary parameters. }
}

\value{
 Returns a list of length one which is the vector of the rates of increase (required by \code{ode}).
}

\references{ 
Scheffer, M., Szabo, S., Gragnani, A., van Nes, E.H., Rinaldi, S.,
Kautsky, N., Norberg, J., Roijackers, R.M.M. and Franken,
R.J.M. (2003) Floating plant dominance as a stable state. \emph{Proceeding of
  the National Academy of Sciences, U.S.A.}, \bold{100}, 4040--4045.

Stevens, M.H.H. (2009) \emph{A Primer of Ecology with R}. Use R! Series. Springer.
 }
\author{ Hank Stevens <HStevens@muohio.edu> }
\seealso{\code{\link{lvcompg}}, \code{\link{igp}}}
\examples{
p <- c(N=2.5, as=0.01, af=0.01, b=0.02, qs=0.075, qf=0.005, 
       hs=0, hf=0.2, ls=0.05, lf=0.05, rs=0.5, rf=0.5, W=0)
t <- 1:200
Initial <- c(F=10, S=10)
S.out1 <- ode(Initial, t, scheffer, p)
matplot(t, S.out1[,-1], type='l')
legend('right', c("F", "S"), lty=1:2, col=1:2, bty='n')
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ methods }