"lmomcoexta.cluster3d" <-
function(lmomdataframe, lcoratdiagram = NULL, site=NA, group=NA,
         file="lmomcoexta.cluster3d.pdf", dopdf=FALSE,
         xlim=c(0,1), ylim=c(-1,1), zlim=c(-1,1), draworigins=TRUE,
         xlab=NULL, ylab=NULL, zlab=NULL,
         arecomoments=TRUE, editlcoratdiagram=TRUE,
         byangle=5, begangle=0, endangle=360, type="p", k=NULL, h=NULL,
         verbose=TRUE, ...) {

   if(is.null(k) && is.null(h)) stop("either 'k' or 'h' must be specified")
   if(is.null(lcoratdiagram)  ) stop("Need lcoratdiagram() object")

   z <- lcoratdiagram
   sites  <- z$sites
   groups <- z$groups

   if(arecomoments) {
     if(is.null(xlab)) xlab <- "L-CORRELATION"
     if(is.null(ylab)) ylab <- "L-COSKEW"
     if(is.null(zlab)) zlab <- "L-COKURTOSIS"
   } else {
     if(is.null(xlab)) xlab <- "L-CV"
     if(is.null(ylab)) ylab <- "L-SKEW"
     if(is.null(zlab)) zlab <- "L-KURTOSIS"
   }

   lmr <- lmomdataframe
   hc  <- hclust(dist(lmr))
   if(! is.null(k)) {
     if(k > 8) stop("k > 8, but only colors 1:8 are implemented")
     col <- cutree(hc, k=1:k)[,k]
   } else {
     col <- cutree(hc, h=h)
     if(max(col) > 8) stop("h is too small, accessing more than 8 colors, but only colors 1:8 are implemented")
   }

   if(verbose) message("lmomcoexta.cluster3d: angle=", appendLF = FALSE)

   if(dopdf) {
     pdf(file)
     for(x in seq(begangle, endangle, by=byangle)) {
        f3d <- scatterplot3d::scatterplot3d(lmr, type=type, xlim=xlim, ylim=ylim, zlim=zlim,
                             color=col, angle=x, pch=16,
                             xlab=xlab, ylab=ylab, zlab=zlab)
        if(draworigins) {
          oc <- rgb(0,0,0,.5)
          t2s <- seq(min(xlim), max(xlim), by=.2)
          for(t2 in t2s[2:length(t2s)]) {
            f3d$points3d(c(t2,t2), ylim, c(0,0), lty=1, col=oc, type="l")
            f3d$points3d(c(t2,t2), c(0,0), zlim, lty=1, col=oc, type="l")
          }
          for(t3 in seq(-.8,.8, by=.2)) {
            f3d$points3d(xlim, c(t3,t3), c(0,0), lty=1, col=oc, type="l")
          }
          for(t4 in seq(-.8,.8, by=.2)) {
            f3d$points3d(xlim, c(0,0), c(t4,t4), lty=1, col=oc, type="l")
          }
          f3d$points3d(xlim, c(-1,-1), c(0,0), lty=1, col=oc, type="l")
          f3d$points3d(xlim, c(1,1),   c(0,0), lty=1, col=oc, type="l")
          f3d$points3d(rep(min(xlim),2), c(-1,1), c(0,0), lty=1, col=oc,  type="l")
          f3d$points3d(rep(max(xlim),2), c(-1,1), c(0,0), lty=1, col=oc,  type="l")
          f3d$points3d(xlim, c(0,0),     c(-1,-1),     lty=1, col=oc, type="l")
          f3d$points3d(xlim, c(0,0),     c(1,1),       lty=1, col=oc, type="l")
          f3d$points3d(rep(min(xlim),2), c(0,0), zlim, lty=1, col=oc, type="l")
          f3d$points3d(rep(max(xlim),2), c(0,0), zlim, lty=1, col=oc, type="l")
        }
        if(! is.na(group)) {
           f3d$points3d(lmr[groups == group,], pch=1, cex=1.5,
                        col=col[groups == group])
        }
        message(x,", ", appendLF = FALSE)
     }
     dev.off()
   }
   if(verbose) message("... done")

   if(editlcoratdiagram) {
      z$color <- col
      return(z)
   } else {
      return(data.frame(site=sites, group=groups, color=col))
   }
}
