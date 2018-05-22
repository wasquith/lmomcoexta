\name{lcoratdiagram}
\alias{lcoratdiagram}
\title{ Create a Special Object of L-moments and L-comoments }
\description{
This function computes an \code{lcoratdiagram} object used by other functions of this package. The object is specially designed for exploratory L-moment analysis of bivariate data. The \code{lcoratdiagram} object is a complex assembly of L-moments, L-comoments, Spearman's Rho, Kendall's Tau, and other elements. The function provides for specialized large-scale visualization of data. The core data are identified by a site variable (\code{site}), the paired x and y values (\code{x}, \code{y}) for each site, and an optional grouping or classification variable (\code{group}). The x and y data could also be either left or right censored, which are identified by indicator variable. A multi-page PDF output file (\code{xyfile}) of the bivariate relation between x and y and corresponding plot of nonexceedance probabilities by plotting position is created. The limits of the bivariate plot (\code{xlim} and \code{ylim}) can be controlled as well as the axis transformation (\code{xtrans}, \code{ytrans}). The labels of the bivariate plot also can be set (\code{xlab}, \code{ylab}). The function computes the first 5 L-moments and L-comoments, therefore dataset having less than 6 values are skipped by default (\code{n.to.next}), but the user is welcome to increase this amount. The function also supports the simulation of bivariate normal distributions by a simulation count (\code{nsim}). These simulated values can be useful for the visualization (\code{\link{plotlcoratdiagram}} or \code{\link{plotlratdiagram}}) of the uncertainty of the relations between L-moment ratios and L-comoment ratios.
}
\usage{
lcoratdiagram(x, y, site, group=NULL,
              xlc=NULL, ylc=NULL, xrc=NULL, yrc=NULL,
              nsim=0, n.to.next=6, verbose=TRUE,
              xyfile="XY.pdf", xlim=NULL, ylim=NULL,
              xtrans=function(x) { return(x) },
              ytrans=function(y) { return(y) },
              xlab = "X VARIABLE",
              ylab = "Y VARIABLE", ...)
}
\arguments{
  \item{lcoratdiagram}{A object from \code{lcoratdiagmra()}.}
  \item{site}{The site id for each data value.}
  \item{group}{An optional grouping value, which can be used to enhance the visualization of the \code{\link{plotlcoratdiagram}} and \code{\link{plotlratdiagram}} options.}
  \item{x}{.}
  \item{y}{.}
  \item{xlc}{.}
  \item{ylc}{.}
  \item{xrc}{Right censoring of x variable by indicator variable (not implemented).}
  \item{yrc}{Right censoring by y variable indicator variable (not implemented).}
  \item{nsim}{.}
  \item{n.to.next}{.}
  \item{verbose}{A logical that controls whether \code{message} updates on status are shown.}
  \item{xyfile}{.}
  \item{xlim}{.}
  \item{ylim}{.}
  \item{xtrans}{.}
  \item{ytrans}{.}
  \item{xlab}{.}
  \item{ylab}{.}
  \item{...}{Additional arguments to pass.}
}
\value{
  An R \code{list} is returned.
}
\author{ W.H. Asquith}

\seealso{ \code{\link{plotlcoratdiagram}}, \code{\link{plotlratdiagram}} }
\examples{
# 
}
\keyword{distribution}
