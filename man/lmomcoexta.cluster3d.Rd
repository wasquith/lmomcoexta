\name{lmomcoexta.cluster3d}
\alias{lmomcoexta.cluster3d}
\title{ Plot 3-D L-moment or L-comoment Ratio Diagrams and Return the Color of Clustering }
\description{
Plot the L-moment Ratio Diagrams from the \code{lcoratdiagram} object from the \code{lcoratdiagram()} function.
}
\usage{
lmomcoexta.cluster3d(lmomdataframe, lcoratdiagram = NULL,
                     site=NA, group=NA,
                     file="lmomcoexta.cluster3d.pdf",
                     dopdf=FALSE,
                     xlim=c(0,1), ylim=c(-1,1), zlim=c(-1,1),
                     draworigins=TRUE,
                     xlab=NULL, ylab=NULL, zlab=NULL,
                     arecomoments=TRUE, editlcoratdiagram=TRUE,
                     byangle=5, begangle=0, endangle=360, type="p",
                     k=NULL, h=NULL, verbose=TRUE, ...)
}
\arguments{
  \item{lmomdataframe}{An x,y,z data frame of L-moments or L-comoments or other for the clustering algorithm.}
  \item{lcoratdiagram}{An object from \code{\link{lcoratdiagram}} function.}
  \item{site}{An optional specific site to identify on the plot. NOT IMPLEMENTED.}
  \item{group}{An optional group of interest to specially plot, this is matched against those groups provide by similarly named argument to \code{\link{lcoratdiagram}}.}
  \item{file}{The PDF file to write in the current working directory \code{getwd()}.}
  \item{dopdf}{Should the PDF file be rendered.}
  \item{xlim}{The limits of the x-axis for the plot. The default implies L-correlation.}
  \item{ylim}{The limits of the y-axis for the plot. The default implies L-skew/L-coskew.}
  \item{zlim}{The limits of the y-axis for the plot. The default implies L-kurtosis/L-cokurtosis.}
  \item{draworigins}{Logical to draw the origin planes on the plot.}
  \item{arecomoments}{A logical on whether the L-comoments are incoming. Controls the axis labels.}
  \item{xlab}{An optional label for the x-axis, defaults are internally applied.}
  \item{ylab}{An optional label for the y-axis, defaults are internally applied.}
  \item{zlab}{An optional label for the z-axis, defaults are internally applied.}
  \item{editlcoratdiagram}{A logical that if true the \code{lcoratdiagram} object is edited to include the \code{color} of the clustering algorithm. If false, a simple data frame of the site, group, and color from the \code{lcoratdiagram} are returned.}
  \item{byangle}{The stepping angle for the repeated rendering of the PDF.}
  \item{begangle}{The beginning angle for the repeated rendering of the PDF.}
  \item{endangle}{The ending angle for the repeated rendering of the PDF.}
  \item{type}{The same \code{type} argument of the \code{scatterplot3d} function from the \pkg{scatterplot3d} package.}
  \item{k}{The same \code{k} argument of the \code{scatterplot3d} function from the \pkg{scatterplot3d} package.}
  \item{h}{The same \code{h} argument of the \code{scatterplot3d} function from the \pkg{scatterplot3d} package.}
  \item{verbose}{A logical that controls whether \code{message} updates on status are shown.}
  \item{...}{Additional arguments not yet used.}
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
