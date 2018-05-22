\name{lcoratdiagram.discordance}
\alias{lcoratdiagram.discordance}
\title{ The Discordance of L-moments and L-comoments }
\description{
Compute the discorance of sample L-moments and L-comoments for a \code{lcoratdiagram} object. 
}
\usage{
lcoratdiagram.discordance(lcoratdiagram=NULL, verbose=TRUE, ...)
}
\arguments{
  \item{lcoratdiagram}{A object from \code{lcoratdiagmra()}.}
  \item{verbose}{A logical that controls whether \code{message} updates on status are shown.}
  \item{...}{Additional arguments to pass to \code{lmrdiagram()} of the \pkg{lmomco} package.}
}
\value{
  An R \code{list} is returned.

  \item{Dx}{Vector of the L-moments. First element is
\eqn{\hat{\lambda}^{(0,0)}_1}, second element is \eqn{\hat{\lambda}^{(0,0)}_2}, and so on.}
  \item{Dy}{Vector of the L-moment ratios. Second element is 
\eqn{\hat{\tau}^{(0,0)}}, third element is \eqn{\hat{\tau}^{(0,0)}_3} and so on.}
  \item{Dxy}{Level of symmetrical trimming used in the computation, which will equal \code{NULL} if asymmetrical trimming was used.}
  \item{Dyx}{Level of left-tail trimming used in the computation.}
}
\note{
This function computes the L-moments through the generalization of the TL-moments (\code{\link{TLmoms}}). In fact, this function calls the default TL-moments with no trimming of the sample. This function is equivalent to \code{\link{lmom.ub}}, but returns a different data structure.
}
\references{
see references in \code{lmrdiscord()} of the \pkg{lmomco} package.
}
\author{ W.H. Asquith}

\seealso{ \code{\link{lcoratdiagram}} }
\examples{
# 
}
\keyword{distribution}
