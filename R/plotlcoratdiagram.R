"plotlcoratdiagram" <-
function(lcoratdiagram = NULL, group=NA,
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
         verbose=TRUE, ...) {

   if(is.null(lcoratdiagram)) stop("Need lcoratdiagram() object")
   if(showhiDsite) showsite <- FALSE

   if(verbose) message("\n\n*** plotlcoratdiagram(), initiation")


   DIS <- lcoratdiagram.discordance(lcoratdiagram=lcoratdiagram, sort=FALSE, ...)

   z <- lcoratdiagram

   W <- z$record.lengths
   TAU <- z$kendall.tau
   RHO <- z$spearman.rho
   sites <- z$sites
   groups <- z$groups

   T2.12 <- z$lcomoments$T2.12
   T2.21 <- z$lcomoments$T2.21
   T3.12 <- z$lcomoments$T3.12
   T3.21 <- z$lcomoments$T3.21
   T4.12 <- z$lcomoments$T4.12
   T4.21 <- z$lcomoments$T4.21

   rT2.12 <- z$regional.means$T2.12
   rT2.21 <- z$regional.means$T2.21
   rT3.12 <- z$regional.means$T3.12
   rT3.21 <- z$regional.means$T3.21
   rT4.12 <- z$regional.means$T4.12
   rT4.21 <- z$regional.means$T4.21

   T2T3cov.12 <- z$covariances$t2t3cov.12
   T2T3cov.21 <- z$covariances$t2t3cov.21
   T3T4cov.12 <- z$covariances$t3t4cov.12
   T3T4cov.21 <- z$covariances$t3t4cov.21
   T2T4cov.12 <- z$covariances$t2t4cov.12
   T2T4cov.21 <- z$covariances$t2t4cov.21

   R <- range(W)
   lwds <- ((W-R[1])/R[2])
   lwds[lwds == 0] <- min(lwds[lwds != 0])

   if(verbose) message(" Working on PDF output of L-comoment Diagrams")

   pdf(file=pdft2t3file, useDingbats=FALSE)
   plot(c(0,0), c(0,0), type="n",
         xlim=xt2t3lim, ylim=yt2t3lim,
         xlab="L-CORRELATION", ylab="L-COSKEW")
   if(showorigin) {
      lines(c(-2, 2), c( 0, 0), lty=2)
      lines(c( 0, 0), c(-2, 2), lty=2)
   }
   if(showbivarsim & ! is.null(z$bivarsim)) {
      T2T3bivn.12  <- z$bivarsim$t2t3bivn.12
      T2T3bivn.21  <- z$bivarsim$t2t3bivn.21
      if(doxy) points(T2T3bivn.12[,1],T2T3bivn.12[,2], pch=16, col=rgb(1,0,0,.4), cex=0.5)
      if(doyx) points(T2T3bivn.21[,1],T2T3bivn.21[,2], pch=16, col=rgb(0,1,0,.4), cex=0.5)
   }

   if(doxy) lines(ellipse::ellipse(T2T3cov.12, centre=c(rT2.12, rT3.12), level=0.50), col=2)
   if(doyx) lines(ellipse::ellipse(T2T3cov.21, centre=c(rT2.21, rT3.21), level=0.50), col=3)
   if(doxy) lines(ellipse::ellipse(T2T3cov.12, centre=c(rT2.12, rT3.12), level=0.90), col=2, lty=2)
   if(doyx) lines(ellipse::ellipse(T2T3cov.21, centre=c(rT2.21, rT3.21), level=0.90), col=3, lty=2)

   for(i in 1:length(W)) {
      if(! is.na(group) & groups[i] == group) {
         if(doxy) points(T2.12[i], T3.12[i], pch=12, cex=1.25, lwd=1.5, col=rgb(0.25,0.25,0.75,1))
         if(doyx) points(T2.21[i], T3.21[i], pch=12, cex=1.25, lwd=1.5, col=rgb(0.25,0.25,0.75,1))
      }
      if(RHO[i] > rhothres) {
         if(doxy) points(T2.12[i], T3.12[i], pch=2, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
         if(doyx) points(T2.21[i], T3.21[i], pch=2, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
      }
      if(TAU[i] > tauthres) {
         if(doxy) points(T2.12[i], T3.12[i], pch=6, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
         if(doyx) points(T2.21[i], T3.21[i], pch=6, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
      }
      if(doxy &DIS$Dxy$isD[i]) points(T2.12[i], T3.12[i], cex=1.5, lwd=.75)
      if(DIS$Dx$isD[i])  {
         if(doxy) points(T2.12[i], T3.12[i], cex=1.5, lwd=.75, pch=3)
         if(doyx) points(T2.21[i], T3.21[i], cex=1.5, lwd=.75, pch=3)
      }
      if(doyx & DIS$Dyx$isD[i]) points(T2.21[i], T3.21[i], cex=1.5, lwd=.75)
      if(DIS$Dy$isD[i])  {
         if(doxy) points(T2.12[i], T3.12[i], cex=1.5, lwd=.75, pch=4)
         if(doyx) points(T2.21[i], T3.21[i], cex=1.5, lwd=.75, pch=4)
      }
      if(barbells) lines(c(T2.12[i], T2.21[i]), c(T3.12[i], T3.21[i]),
                         lwd=barbellwd*lwds[i], col=rgb(0,0,0,lwds[i]^(barbellpw)), pch=16)
      if(doxy) points(T2.12[i], T3.12[i], pch=21, bg=rgb(1,0,0,.8), cex=.75)
      if(doyx) points(T2.21[i], T3.21[i], pch=21, bg=rgb(0,1,0,.8), cex=.75)
      if(showsite) {
         if(doxy) text(T2.12[i], T3.12[i], labels=sites[i], cex=0.5, adj=c(0,-1))
         if(doyx) text(T2.21[i], T3.21[i], labels=sites[i], cex=0.5, adj=c(0, 1))
      }
      else if(showhiDsite) {
         if(DIS$Dxy$isD[i] || DIS$Dx$isD[i] || DIS$Dyx$isD[i] || DIS$Dy$isD[i]) {
            if(doxy) text(T2.12[i], T3.12[i], labels=sites[i], cex=0.5, adj=c(0,-1))
            if(doyx) text(T2.21[i], T3.21[i], labels=sites[i], cex=0.5, adj=c(0, 1))
         }
      }
   }
   if(doxy & doyx) lines(c(rT2.12, rT2.21), c(rT3.12,rT3.21), col=2)
   if(doxy) points(rT2.12, rT3.12, pch=21, bg=rgb(1,0,0,1), cex=2)
   if(doyx) points(rT2.21, rT3.21, pch=21, bg=rgb(0,1,0,1), cex=2)
   dev.off()
   if(verbose) message(" PDF coT2T3 file=",pdft2t3file," in ",getwd())

   pdf(file=pdft3t4file)
   plot(c(0,0), c(0,0), type="n",
         xlim=xt3t4lim, ylim=yt3t4lim,
         xlab="L-COSKEW", ylab="L-COKURTOSIS")
   if(showorigin) {
      lines(c(-2, 2), c( 0, 0), lty=2)
      lines(c( 0, 0), c(-2, 2), lty=2)
   }
   if(showbivarsim & ! is.null(z$bivarsim)) {
      T3T4bivn.12  <- z$bivarsim$t3t4bivn.12
      T3T4bivn.21  <- z$bivarsim$t3t4bivn.21
      if(doxy) points(T3T4bivn.12[,1],T3T4bivn.12[,2], pch=16, col=rgb(1,0,0,.4), cex=0.5)
      if(doyx) points(T3T4bivn.21[,1],T3T4bivn.21[,2], pch=16, col=rgb(0,1,0,.4), cex=0.5)
   }

   if(doxy) lines(ellipse::ellipse(T3T4cov.12, centre=c(rT3.12, rT4.12), level=0.50), col=2)
   if(doyx) lines(ellipse::ellipse(T3T4cov.21, centre=c(rT3.21, rT4.21), level=0.50), col=3)
   if(doxy) lines(ellipse::ellipse(T3T4cov.12, centre=c(rT3.12, rT4.12), level=0.90), col=2, lty=2)
   if(doyx) lines(ellipse::ellipse(T3T4cov.21, centre=c(rT3.21, rT4.21), level=0.90), col=3, lty=2)

   for(i in 1:length(W)) {
      if(! is.na(group) & groups[i] == group) {
         if(doxy) points(T3.12[i], T4.12[i], pch=12, cex=1.25, lwd=1.5, col=rgb(0.25,0.25,0.75,1))
         if(doyx) points(T3.21[i], T4.21[i], pch=12, cex=1.25, lwd=1.5, col=rgb(0.25,0.25,0.75,1))
      }
      if(RHO[i] > rhothres) {
         if(doxy) points(T3.12[i], T4.12[i], pch=2, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
         if(doyx) points(T3.21[i], T4.21[i], pch=2, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
      }
      if(TAU[i] > tauthres) {
         if(doxy) points(T3.12[i], T4.12[i], pch=6, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
         if(doyx) points(T3.21[i], T4.21[i], pch=6, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
      }
      if(doxy & DIS$Dxy$isD[i]) points(T3.12[i], T4.12[i], cex=1.5, lwd=.75)
      if(DIS$Dx$isD[i])  {
         if(doxy) points(T3.12[i], T4.12[i], cex=1.5, lwd=.75, pch=3)
         if(doyx) points(T3.21[i], T4.21[i], cex=1.5, lwd=.75, pch=3)
      }

      if(doyx & DIS$Dyx$isD[i]) points(T3.21[i], T4.21[i], cex=1.5, lwd=.75)
      if(DIS$Dy$isD[i])  {
         if(doxy) points(T3.12[i], T4.12[i], cex=1.5, lwd=.75, pch=4)
         if(doyx) points(T3.21[i], T4.21[i], cex=1.5, lwd=.75, pch=4)
      }
      if(barbells) lines(c(T3.12[i], T3.21[i]), c(T4.12[i], T4.21[i]),
                         lwd=barbellwd*lwds[i], col=rgb(0,0,0,lwds[i]^(barbellpw)), pch=16)
      if(doxy) points(T3.12[i], T4.12[i], pch=21, bg=rgb(1,0,0,.8), cex=.75)
      if(doyx) points(T3.21[i], T4.21[i], pch=21, bg=rgb(0,1,0,.8), cex=.75)
      if(showsite) {
         if(doxy) text(T3.12[i], T4.12[i], labels=sites[i], cex=0.5, adj=c(0,-1))
         if(doyx) text(T3.21[i], T4.21[i], labels=sites[i], cex=0.5, adj=c(0, 1))
      } else if(showhiDsite) {
         if(DIS$Dxy$isD[i] || DIS$Dx$isD[i] || DIS$Dyx$isD[i] || DIS$Dy$isD[i]) {
            if(doxy) text(T3.12[i], T4.12[i], labels=sites[i], cex=0.5, adj=c(0,-1))
            if(doyx) text(T3.21[i], T4.21[i], labels=sites[i], cex=0.5, adj=c(0, 1))
         }
      }
   }
   if(doxy & doyx) lines(c(rT3.12, rT3.21), c(rT4.12,rT4.21), col=2)
   if(doxy) points(rT3.12, rT4.12, pch=21, bg=rgb(1,0,0,1), cex=2)
   if(doyx) points(rT3.21, rT4.21, pch=21, bg=rgb(0,1,0,1), cex=2)
   dev.off()
   if(verbose) message(" PDF coT3T4 file=",pdft3t4file," in ",getwd())


   pdf(file=pdft2t4file)
   plot(c(0,0), c(0,0), type="n",
         xlim=xt2t4lim, ylim=yt2t4lim,
         xlab="L-CORRELATION", ylab="L-COKURTOSIS")
   if(showorigin) {
      lines(c(-2, 2), c( 0, 0), lty=2)
      lines(c( 0, 0), c(-2, 2), lty=2)
   }
   if(showbivarsim & ! is.null(z$bivarsim)) {
      T2T4bivn.12  <- z$bivarsim$t2t4bivn.12
      T2T4bivn.21  <- z$bivarsim$t2t4bivn.21
      if(doxy) points(T2T4bivn.12[,1],T2T4bivn.12[,2], pch=16, col=rgb(1,0,0,.4), cex=0.5)
      if(doyx) points(T2T4bivn.21[,1],T2T4bivn.21[,2], pch=16, col=rgb(0,1,0,.4), cex=0.5)
   }

   if(doxy) lines(ellipse::ellipse(T2T4cov.12, centre=c(rT2.12, rT4.12), level=0.50), col=2)
   if(doyx) lines(ellipse::ellipse(T2T4cov.21, centre=c(rT2.21, rT4.21), level=0.50), col=3)
   if(doxy) lines(ellipse::ellipse(T2T4cov.12, centre=c(rT2.12, rT4.12), level=0.90), col=2, lty=2)
   if(doyx) lines(ellipse::ellipse(T2T4cov.21, centre=c(rT2.21, rT4.21), level=0.90), col=3, lty=2)

   for(i in 1:length(W)) {
      if(! is.na(group) & groups[i] == group) {
         if(doxy) points(T2.12[i], T4.12[i], pch=12, cex=1.25, lwd=1.5, col=rgb(0.25,0.25,0.75,1))
         if(doyx) points(T2.21[i], T4.21[i], pch=12, cex=1.25, lwd=1.5, col=rgb(0.25,0.25,0.75,1))
      }
      if(RHO[i] > rhothres) {
         if(doxy) points(T2.12[i], T4.12[i], pch=2, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
         if(doyx) points(T2.21[i], T4.21[i], pch=2, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
      }
      if(TAU[i] > tauthres) {
         if(doxy) points(T2.12[i], T4.12[i], pch=6, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
         if(doyx) points(T2.21[i], T4.21[i], pch=6, cex=1.25, lwd=0.75, col=rgb(0,.2,.8,.7))
      }
      if(doxy & DIS$Dxy$isD[i]) points(T2.12[i], T4.12[i], cex=1.5, lwd=.75)
      if(DIS$Dx$isD[i])  {
         if(doxy) points(T2.12[i], T4.12[i], cex=1.5, lwd=.75, pch=3)
         if(doyx) points(T2.21[i], T4.21[i], cex=1.5, lwd=.75, pch=3)
      }

      if(doyx & DIS$Dyx$isD[i]) points(T2.21[i], T4.21[i], cex=1.5, lwd=.75)
      if(DIS$Dy$isD[i]) {
         if(doxy) points(T2.12[i], T4.12[i], cex=1.5, lwd=.75, pch=4)
         if(doyx) points(T2.21[i], T4.21[i], cex=1.5, lwd=.75, pch=4)
      }

      if(barbells) lines(c(T2.12[i], T2.21[i]), c(T4.12[i], T4.21[i]),
                         lwd=barbellwd*lwds[i], col=rgb(0,0,0,lwds[i]^(barbellpw)), pch=16)
      if(doxy) points(T2.12[i], T4.12[i], pch=21, bg=rgb(1,0,0,.8), cex=.75)
      if(doyx) points(T2.21[i], T4.21[i], pch=21, bg=rgb(0,1,0,.8), cex=.75)
      if(showsite) {
         if(doxy) text(T2.12[i], T4.12[i], labels=sites[i], cex=0.5, adj=c(0,-1))
         if(doyx) text(T2.21[i], T4.21[i], labels=sites[i], cex=0.5, adj=c(0, 1))
      } else if(showhiDsite) {
         if(DIS$Dxy$isD[i] || DIS$Dx$isD[i] || DIS$Dyx$isD[i] || DIS$Dy$isD[i]) {
            if(doxy) text(T2.12[i], T4.12[i], labels=sites[i], cex=0.5, adj=c(0,-1))
            if(doyx) text(T2.21[i], T4.21[i], labels=sites[i], cex=0.5, adj=c(0, 1))
         }
      }
   }
   if(doxy & doyx) lines(c(rT2.12, rT2.21), c(rT4.12,rT4.21), col=2)
   if(doxy) points(rT2.12, rT4.12, pch=21, bg=rgb(1,0,0,1), cex=2)
   if(doyx) points(rT2.21, rT4.21, pch=21, bg=rgb(0,1,0,1), cex=2)
   dev.off()
   if(verbose) message(" PDF coT2T4 file=",pdft2t4file," in ",getwd())

   if(verbose) message("*** plotlcoratdiagram(), done")
}

