"lcoratdiagram.discordance" <-
function(lcoratdiagram=NULL, verbose=TRUE, ...) {

   if(is.null(lcoratdiagram)) stop("Need lcoratdiagram() object")

   z <- lcoratdiagram

   if(verbose) message("\n\n*** lcoratdiagram.discordance(), initiation")

   sites <- z$sites

   T2.12 <- z$lcomoments$T2.12
   T2.21 <- z$lcomoments$T2.21
   T3.12 <- z$lcomoments$T3.12
   T3.21 <- z$lcomoments$T3.21
   T4.12 <- z$lcomoments$T4.12
   T4.21 <- z$lcomoments$T4.21

   xL1 <- z$lmoments$x$L1
   xT2 <- z$lmoments$x$T2
   xT3 <- z$lmoments$x$T3
   xT4 <- z$lmoments$x$T4

   yL1 <- z$lmoments$y$L1
   yT2 <- z$lmoments$y$T2
   yT3 <- z$lmoments$y$T3
   yT4 <- z$lmoments$y$T4

   if(verbose) message("  Discordance table for X variable")
   Dx <- lmomco::lmrdiscord(site=sites, t2=xT2, t3=xT3, t4=xT4, ...)

   if(verbose) message("  Discordance table for Y variable")
   Dy <- lmomco::lmrdiscord(site=sites, t2=yT2, t3=yT3, t4=yT4, ...)

   if(verbose) message("  Discordance table for X wrt Y variable")
   Dxy <- lmomco::lmrdiscord(site=sites, t2=T2.12, t3=T3.12, t4=T4.12, ...)

   if(verbose) message("  Discordance table for Y wrt X variable")
   Dyx <- lmomco::lmrdiscord(site=sites, t2=T2.21, t3=T3.21, t4=T4.21, ...)

   if(verbose) message("*** lcoratdiagram.discordance(), done")

   w <- list(Dx=Dx, Dy=Dy, Dxy=Dxy, Dyx=Dyx)
   return(w)
}


