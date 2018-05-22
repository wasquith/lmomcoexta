\name{plotlcoratdiagram}
\alias{plotlcoratdiagram}
\title{ Plot L-comoment Ratio Diagrams }
\description{
Plot the L-comoment Ratio Diagrams from the \code{lcoratdiagram} object from the \code{lcoratdiagram()} function.
}
\usage{
plotlcoratdiagram(lcoratdiagram = NULL, group=NA,
                  xt2t3lim=c(-1,1), yt2t3lim=c(-1,1),
                  xt3t4lim=c(-1,1), yt3t4lim=c(-1,1),
                  xt2t4lim=c(-1,1), yt2t4lim=c(-1,1),
                  showbivarsim=FALSE, showsite=TRUE, showhiDsite=TRUE,
                  rhothres=Inf, tauthres=Inf,
                  barbells=TRUE, barbellwd=1.2, barbellpw=1/3,
                  showorigin=TRUE,
                  doxy=TRUE, doyx=TRUE,
                  pdft2t3file="lcoratT2T3.pdf",
                  pdft3t4file="lcoratT3T4.pdf",
                  pdft2t4file="lcoratT2T4.pdf",
                  verbose=TRUE, ...)
}
\arguments{
  \item{lcoratdiagram}{An object from \code{\link{lcoratdiagram}} function.}
  \item{group}{An optional group of interest to specially plot, this is matched against those groups provide by similarly named argument to \code{\link{lcoratdiagram}}.}
  \item{xt2t3lim}{.}
  \item{yt2t3lim}{.}
  \item{xt3t4lim}{.}
  \item{yt3t4lim}{.}
  \item{xt2t4lim}{.}
  \item{yt2t4lim}{.}
  \item{showbivarsim}{Show bivariate simulations, which were specified by the \code{nsim} argument of \code{\link{lcoratdiagram}}.}
  \item{showsite}{A logical on whether site ids are shown in the plots.}
  \item{showhiDsite}{A logical on whether ids for highly discordant sites are shown. If \code{TRUE}, then \code{showsite} is taken as \code{FALSE}.}
  \item{rhothres}{A threshold on Spearman's Rho. For Rho values greater than this threshold, special plotting operations are made.}
  \item{tauthres}{A threshold on Kendall's Tau. For Tau values greater than this threshold, special plotting operations are made.}
  \item{barbells}{A logical toggling whether barbells are drawn.}
  \item{barbellwd}{A factor multiplier on the barbell thicknesses, which are also dependent on the number of samples per site.}
  \item{barbellpw}{An exponent on the alpha of the \code{rgb} function for the color. This exponent should be less than 1 and increasingly less than 1 darkens the barbells.}  
  \item{showorigin}{A logical toggling the drawing of the origin lines for the respective plots.}
  \item{doxy}{Show only symbols for the x with respect to y variable.}
  \item{doyx}{Show only symbols for the y with respect to x variable.}
  \item{pdft2t3file}{.}
  \item{pdft3t4file}{.}
  \item{pdft2t4file}{.}
  \item{verbose}{A logical that controls whether \code{message} updates on status are shown.}
  \item{...}{Additional arguments to pass to the \code{\link{lcoratdiagram.discordance}} function.}
}
\value{
  Stuff happens.
}
\author{ W.H. Asquith}
\seealso{ \code{\link{lcoratdiagram}} }
\examples{
# 
}
\keyword{distribution}
