"plotlratdiagram" <-
function(lcoratdiagram = NULL, group=NA,
         xt2t3lim=c(0,1),  yt2t3lim=c(-1,1),
         xt3t4lim=c(-1,1), yt3t4lim=c(-1,1),
         xt2t4lim=c(0,1),  yt2t4lim=c(-1,1),
         showbivarsim=FALSE, showsite=TRUE, showhiDsite=TRUE,
         rhothres=Inf, tauthres=Inf,
         barbells=TRUE, barbellwd=1.2, barbellpw=1/3,
         showorigin=TRUE,
         dox=TRUE, doy=TRUE,
         pdft2t3file="lratT2T3.pdf",
         pdft3t4file="lratT3T4.pdf",
         pdft2t4file="lratT2T4.pdf",
         verbose=TRUE, ...) {

   if(is.null(lcoratdiagram)) stop("Need lcoratdiagram() object")
   if(showhiDsite) showsite <- FALSE

   if(verbose) message("\n\n*** plotlratdiagram(), initiation")

   DIS <- lcoratdiagram.discordance(lcoratdiagram=lcoratdiagram, sort=FALSE, ...)

   z <- lcoratdiagram

   W <- z$record.lengths
   TAU <- z$kendall.tau
   RHO <- z$spearman.rho
   sites <- z$sites
   groups <- z$groups

   x <- z$lmoments$x
   xL1 <- x$L1
   xT2 <- x$T2
   xT3 <- x$T3
   xT4 <- x$T4

   y <- z$lmoments$y
   yL1 <- y$L1
   yT2 <- y$T2
   yT3 <- y$T3
   yT4 <- y$T4

   x <- z$regional.means$x
   rxL1 <- x$L1
   rxT2 <- x$T2
   rxT3 <- x$T3
   rxT4 <- x$T4

   y <- z$regional.means$y
   ryL1 <- y$L1
   ryT2 <- y$T2
   ryT3 <- y$T3
   ryT4 <- y$T4

   xT2T3cov <- z$covariances$x$t2t3cov;
   xT3T4cov <- z$covariances$x$t3t4cov;
   xT2T4cov <- z$covariances$x$t2t4cov;

   yT2T3cov <- z$covariances$y$t2t3cov;
   yT3T4cov <- z$covariances$y$t3t4cov;
   yT2T4cov <- z$covariances$y$t2t4cov;


   xT2T3bivn <- z$bivarsim$x$t2t3
   xT3T4bivn <- z$bivarsim$x$t3t4
   xT2T4bivn <- z$bivarsim$x$t2t4

   yT2T3bivn <- z$bivarsim$y$t2t3
   yT3T4bivn <- z$bivarsim$y$t3t4
   yT2T4bivn <- z$bivarsim$y$t2t4


   R <- range(W)
   lwds <- ((W-R[1])/R[2])
   lwds[lwds == 0] <- min(lwds[lwds != 0])

   if(verbose) message(" Working on PDF output of L-comoment Diagrams")

   pdf(file=pdft2t3file)
   plot(c(0,0), c(0,0), type="n",
         xlim=xt2t3lim, ylim=yt2t3lim,
         xlab="L-CV", ylab="L-SKEW")
   if(showorigin) {
      lines(c(0,1e6), c(0,0), lty=2)
   }
   if(showbivarsim & ! is.null(z$bivarsim)) {
      if(dox) points(xT2T3bivn[,1],xT2T3bivn[,2], pch=16, col=rgb(1,0,0,.4), cex=0.5)
      if(doy) points(yT2T3bivn[,1],yT2T3bivn[,2], pch=16, col=rgb(0,1,0,.4), cex=0.5)
   }

   if(dox) lines(ellipse::ellipse(xT2T3cov, centre=c(rxT2, rxT3), level=0.50), col=2)
   if(doy) lines(ellipse::ellipse(yT2T3cov, centre=c(ryT2, ryT3), level=0.50), col=3)
   if(dox) lines(ellipse::ellipse(xT2T3cov, centre=c(rxT2, rxT3), level=0.90), col=2, lty=2)
   if(doy) lines(ellipse::ellipse(yT2T3cov, centre=c(ryT2, ryT3), level=0.90), col=3, lty=2)

   for(i in 1:length(W)) {
      if(! is.na(group) & groups[i] == group) {
         if(dox) points(xT2[i], xT3[i], pch=12, cex=1.25, lwd=1.5, col=rgb(0.25,0.25,0.75,1))
         if(doy) points(yT2[i], yT3[i], pch=12, cex=1.25, lwd=1.5, col=rgb(0.25,0.25,0.75,1))
      }
      if(RHO[i] > rhothres) {
         if(dox) points(xT2[i], xT3[i], pch=2, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
         if(doy) points(yT2[i], yT3[i], pch=2, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
      }
      if(TAU[i] > tauthres) {
         if(dox) points(xT2[i], xT3[i], pch=6, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
         if(doy) points(yT2[i], yT3[i], pch=6, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
      }
      if(DIS$Dxy$isD[i]) {
         if(dox) points(xT2[i], xT3[i], cex=1.5, lwd=.75, pch=3)
         if(doy) points(yT2[i], yT3[i], cex=1.5, lwd=.75, pch=3)
      }
      if(DIS$Dx$isD[i])  {
         if(dox) points(xT2[i], xT3[i], cex=1.5, lwd=.75)
      }
      if(DIS$Dyx$isD[i]) {
         if(dox) points(xT2[i], xT3[i], cex=1.5, lwd=.75, pch=4)
         if(doy) points(yT2[i], yT3[i], cex=1.5, lwd=.75, pch=4)
      }
      if(DIS$Dy$isD[i])  {
         if(doy) points(yT2[i], yT3[i], cex=1.5, lwd=.75)
      }
      if(barbells) lines(c(xT2[i], yT2[i]), c(xT3[i], yT3[i]),
                         lwd=barbellwd*lwds[i], col=rgb(0,0,0,lwds[i]^(barbellpw)), pch=16)
      if(dox) points(xT2[i], xT3[i], pch=21, bg=rgb(1,0,0,.8), cex=.75)
      if(doy) points(yT2[i], yT3[i], pch=21, bg=rgb(0,1,0,.8), cex=.75)
      if(showsite) {
         if(dox) text(xT2[i], xT3[i], labels=sites[i], cex=0.5, adj=c(0,-1))
         if(doy) text(yT2[i], yT3[i], labels=sites[i], cex=0.5, adj=c(0, 1))
      } else if(showhiDsite) {
         if(DIS$Dxy$isD[i] | DIS$Dx$isD[i] | DIS$Dyx$isD[i] | DIS$Dy$isD[i]) {
            if(dox) text(xT2[i], xT3[i], labels=sites[i], cex=0.5, adj=c(0,-1))
            if(doy) text(yT2[i], yT3[i], labels=sites[i], cex=0.5, adj=c(0, 1))
         }
      }
   }
   if(dox & doy) lines(c(rxT2, ryT2), c(rxT3,ryT3), col=2)
   if(dox) points(rxT2, rxT3, pch=21, bg=rgb(1,0,0,1), cex=2)
   if(doy) points(ryT2, ryT3, pch=21, bg=rgb(0,1,0,1), cex=2)
   dev.off()
   if(verbose) message(" PDF T2T3 file=",pdft2t3file," in ",getwd())

   pdf(file=pdft3t4file)
   plot(c(0,0), c(0,0), type="n",
         xlim=xt3t4lim, ylim=yt3t4lim,
         xlab="L-SKEW", ylab="L-KURTOSIS")
   if(showorigin) {
      lines(c(-2,2), c(0,0), lty=2)
      lines(c(0,0), c(-2,2), lty=2)
   }
   if(showbivarsim & ! is.null(z$bivarsim)) {
      if(dox) points(xT3T4bivn[,1],xT3T4bivn[,2], pch=16, col=rgb(1,0,0,.4), cex=0.5)
      if(doy) points(yT3T4bivn[,1],yT3T4bivn[,2], pch=16, col=rgb(0,1,0,.4), cex=0.5)
   }

   if(dox) lines(ellipse::ellipse(xT3T4cov, centre=c(rxT3, rxT4), level=0.50), col=2)
   if(doy) lines(ellipse::ellipse(yT3T4cov, centre=c(ryT3, ryT4), level=0.50), col=3)
   if(dox) lines(ellipse::ellipse(xT3T4cov, centre=c(rxT3, rxT4), level=0.90), col=2, lty=2)
   if(doy) lines(ellipse::ellipse(yT3T4cov, centre=c(ryT3, ryT4), level=0.90), col=3, lty=2)

   for(i in 1:length(W)) {
      if(! is.na(group) & groups[i] == group) {
         if(dox) points(xT3[i], xT4[i], pch=12, cex=1.25, lwd=1.5, col=rgb(0.25,0.25,0.75,1))
         if(doy) points(yT3[i], yT4[i], pch=12, cex=1.25, lwd=1.5, col=rgb(0.25,0.25,0.75,1))
      }
      if(RHO[i] > rhothres) {
         if(dox) points(xT3[i], xT4[i], pch=2, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
         if(doy) points(yT3[i], yT4[i], pch=2, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
      }
      if(TAU[i] > tauthres) {
         if(dox) points(xT3[i], xT4[i], pch=6, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
         if(doy) points(yT3[i], yT4[i], pch=6, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
      }
      if(DIS$Dxy$isD[i]) {
         if(dox) points(xT3[i], xT4[i], cex=1.5, lwd=.75, pch=3)
         if(doy) points(yT3[i], yT4[i], cex=1.5, lwd=.75, pch=3)
      }
      if(DIS$Dx$isD[i])  {
         if(dox) points(xT3[i], xT4[i], cex=1.5, lwd=.75)
      }
      if(DIS$Dyx$isD[i]) {
         if(dox) points(xT3[i], xT4[i], cex=1.5, lwd=.75, pch=4)
         if(doy) points(yT3[i], yT4[i], cex=1.5, lwd=.75, pch=4)
      }
      if(DIS$Dy$isD[i])  {
         if(doy) points(yT3[i], yT4[i], cex=1.5, lwd=.75)
      }
      if(barbells) lines(c(xT3[i], yT3[i]), c(xT4[i], yT4[i]),
                         lwd=barbellwd*lwds[i], col=rgb(0,0,0,lwds[i]^(barbellpw)), pch=16)
      if(dox) points(xT3[i], xT4[i], pch=21, bg=rgb(1,0,0,.8), cex=.75)
      if(doy) points(yT3[i], yT4[i], pch=21, bg=rgb(0,1,0,.8), cex=.75)
      if(showsite) {
         if(dox) text(xT3[i], xT4[i], labels=sites[i], cex=0.5, adj=c(0,-1))
         if(doy) text(yT3[i], yT4[i], labels=sites[i], cex=0.5, adj=c(0, 1))
      } else if(showhiDsite) {
         if(DIS$Dxy$isD[i] | DIS$Dx$isD[i] | DIS$Dyx$isD[i] | DIS$Dy$isD[i]) {
            if(dox) text(xT3[i], xT4[i], labels=sites[i], cex=0.5, adj=c(0,-1))
            if(doy) text(yT3[i], yT4[i], labels=sites[i], cex=0.5, adj=c(0, 1))
         }
      }
   }
   if(dox & doy) lines(c(rxT3, ryT3), c(rxT4,ryT4), col=2)
   if(dox) points(rxT3, rxT4, pch=21, bg=rgb(1,0,0,1), cex=2)
   if(doy) points(ryT3, ryT4, pch=21, bg=rgb(0,1,0,1), cex=2)
   dev.off()
   if(verbose) message(" PDF T3T4 file=",pdft3t4file," in ",getwd())


   pdf(file=pdft2t4file)
   plot(c(0,0), c(0,0), type="n",
        xlim=xt2t4lim, ylim=yt2t4lim,
        xlab="L-CV", ylab="L-KURTOSIS")
   if(showorigin) {
      lines(c(0,1e6), c(0,0), lty=2)
   }
   if(showbivarsim & ! is.null(z$bivarsim)) {
      if(dox) points(xT2T4bivn[,1],xT2T4bivn[,2], pch=16, col=rgb(1,0,0,.4), cex=0.5)
      if(doy) points(yT2T4bivn[,1],yT2T4bivn[,2], pch=16, col=rgb(0,1,0,.4), cex=0.5)
   }

   if(dox) lines(ellipse::ellipse(xT2T4cov, centre=c(rxT2, rxT4), level=0.50), col=2)
   if(doy) lines(ellipse::ellipse(yT2T4cov, centre=c(ryT2, ryT4), level=0.50), col=3)
   if(dox) lines(ellipse::ellipse(xT2T4cov, centre=c(rxT2, rxT4), level=0.90), col=2, lty=2)
   if(doy) lines(ellipse::ellipse(yT2T4cov, centre=c(ryT2, ryT4), level=0.90), col=3, lty=2)

   for(i in 1:length(W)) {
      if(! is.na(group) & groups[i] == group) {
         if(dox) points(xT2[i], xT4[i], pch=12, cex=1.25, lwd=1.5, col=rgb(0.25,0.25,0.75,1))
         if(doy) points(yT2[i], yT4[i], pch=12, cex=1.25, lwd=1.5, col=rgb(0.25,0.25,0.75,1))
      }
      if(RHO[i] > rhothres) {
         if(dox) points(xT2[i], xT4[i], pch=2, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
         if(doy) points(yT2[i], yT4[i], pch=2, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
      }
      if(TAU[i] > tauthres) {
         if(dox) points(xT2[i], xT4[i], pch=6, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
         if(doy) points(yT2[i], yT4[i], pch=6, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
      }
      if(DIS$Dxy$isD[i]) {
         if(dox) points(xT2[i], xT4[i], cex=1.5, lwd=.75, pch=3)
         if(doy) points(yT2[i], yT4[i], cex=1.5, lwd=.75, pch=3)
      }
      if(DIS$Dx$isD[i])  {
         if(dox) points(xT2[i], xT4[i], cex=1.5, lwd=.75)
      }
      if(DIS$Dyx$isD[i]) {
         if(dox) points(xT2[i], xT4[i], cex=1.5, lwd=.75, pch=4)
         if(doy) points(yT2[i], yT4[i], cex=1.5, lwd=.75, pch=4)
      }
      if(DIS$Dy$isD[i])  {
         if(doy) points(yT2[i], yT4[i], cex=1.5, lwd=.75)
      }
      if(barbells) lines(c(xT2[i], yT2[i]), c(xT4[i], yT4[i]),
                         lwd=barbellwd*lwds[i], col=rgb(0,0,0,lwds[i]^(barbellpw)), pch=16)
      if(dox) points(xT2[i], xT4[i], pch=21, bg=rgb(1,0,0,.8), cex=.75)
      if(doy) points(yT2[i], yT4[i], pch=21, bg=rgb(0,1,0,.8), cex=.75)
      if(showsite) {
         if(dox) text(xT2[i], xT4[i], labels=sites[i], cex=0.5, adj=c(0,-1))
         if(doy) text(yT2[i], yT4[i], labels=sites[i], cex=0.5, adj=c(0, 1))
      } else if(showhiDsite) {
         if(DIS$Dxy$isD[i] | DIS$Dx$isD[i] | DIS$Dyx$isD[i] | DIS$Dy$isD[i]) {
            if(dox) text(xT2[i], xT4[i], labels=sites[i], cex=0.5, adj=c(0,-1))
            if(doy) text(yT2[i], yT4[i], labels=sites[i], cex=0.5, adj=c(0, 1))
         }
      }
   }
   if(dox & doy) lines(c(rxT2, ryT2), c(rxT4,ryT4), col=2)
   if(dox) points(rxT2, rxT4, pch=21, bg=rgb(1,0,0,1), cex=2)
   if(doy) points(ryT2, ryT4, pch=21, bg=rgb(0,1,0,1), cex=2)
   dev.off()
   if(verbose) message(" PDF T2T4 file=",pdft2t4file," in ",getwd())

   if(verbose) message("*** plotlratdiagram(), done")
}

