"lcoratdiagram" <-
function(x, y, site, group=NULL,
         xlc=NULL, ylc=NULL, xrc=NULL, yrc=NULL,
         nsim=0, n.to.next=6, verbose=TRUE,
         xyfile="XY.pdf", xlim=NULL, ylim=NULL,
         xtrans=function(x) { return(x) },
         ytrans=function(y) { return(y) },
         xlab = "X VARIABLE",
         ylab = "Y VARIABLE", ...) {

   if(verbose) message("\n\n*** lcoratdiagram(), initiation")

   z <- list()
   W <- TAU <- RHO <- vector(mode="numeric")
   xL1 <- xT2 <- xT3 <- xT4 <- xT5 <- vector(mode="numeric")
   yL1 <- yT2 <- yT3 <- yT4 <- yT5 <- vector(mode="numeric")
   T2.12 <- T2.21 <- T3.12 <- T3.21 <- T4.12 <- T4.21 <- each.site <- each.group <- W

   if(length(x) != length(y)) stop("x and y must be equal length")
   if(! is.null(group) & length(group) != length(x)) stop("group must be equal in length to x and y")

   if(is.null(xlc)) xlc <- rep(FALSE, length(x))
   if(is.null(ylc)) ylc <- rep(FALSE, length(y))

   pdf(xyfile, useDingbats=FALSE)

   sites <- levels(as.factor(site))
   n.sites <- length(sites)
   i <- j <- 0
   for(my.site in sites) {
      j <- j + 1
      if(is.null(group)) {
        my.group <- NA
      } else {
        my.group <- unique(group[site == my.site])
        if(length(my.group) > 1) {
           stop("HELP")
        }
        my.group <- my.group[1]
      }
      my.y  <-   y[my.site == site]
      my.yc <- ylc[my.site == site]
      my.x  <-   x[my.site == site]
      my.xc <- xlc[my.site == site]
      my.n  <- length(my.x)

      if(my.n < n.to.next) {
        if(verbose) message("Skipping site j=",j," or site=",my.site,": to few data points (<",n.to.next,")")
        next;
      }

      plotting.x <- xtrans(my.x)
      plotting.y <- ytrans(my.y)

      ifelse(is.null(xlim), my.xlim <- range(plotting.x), my.xlim <- xlim)
      ifelse(is.null(ylim), my.ylim <- range(plotting.y), my.ylim <- ylim)

      i <- i + 1

      if(verbose) message("  Processing site i=",i," or site=",my.site," a member of group=",my.group," of total=",n.sites," sites")
      each.site[i] <- my.site
      each.group[i] <- my.group

      W[i]   <- my.n
      tau <- cor(my.x, my.y, method="kendall")
      rho <- cor(my.x, my.y, method="spearman")
      TAU[i] <- tau
      RHO[i] <- rho
      if(verbose) message("    Kendall's Tau and Spearman's Rho, done")

      layout(matrix(1:2, nrow=2))
      plot(plotting.x,plotting.y,
           xlim=my.xlim, ylim=my.ylim,
           xlab=xlab, ylab=ylab)
      mtext(paste(c(my.site,"   RHO=",round(rho,digits=3),"; TAU=",round(tau,digits=3)),sep="", collapse=""))

      plot(lmomco::pp(my.x, sort=FALSE), lmomco::pp(my.y, sort=FALSE),
           xlab="X NONEXCEEDANCE PROBABILITY",
           ylab="Y NONEXCEEDANCE PROBABILITY",
           xlim=c(0,1), ylim=c(0,1))
      mtext(paste(c(my.site,"   RHO=",round(rho,digits=3),"; TAU=",round(tau,digits=3)),sep="", collapse=""))


      my.y.lmr <- lmomco::fliplmoms(lmomco::lmomsRCmark(my.y, rcmark=my.yc, flip=TRUE))
      my.x.lmr <- lmomco::fliplmoms(lmomco::lmomsRCmark(my.x, rcmark=my.xc, flip=TRUE))
      if(verbose) message("    Flipped left-censored L-moments, done")


      lcom <- NA
      lcom <- lmomco::lcomoms2(data.frame(X=my.x, Y=my.y), nmom=4)
      mT2.12 <- lcom$T2[1,2]
      mT2.21 <- lcom$T2[2,1]
      mT3.12 <- lcom$T3[1,2]
      mT3.21 <- lcom$T3[2,1]
      mT4.12 <- lcom$T4[1,2]
      mT4.21 <- lcom$T4[2,1]
      if(verbose) message("    L-comoments (r=2,3,4), done")

      xL1[i] <- my.x.lmr$lambdas[1]
      xT2[i] <- my.x.lmr$lambdas[2]/my.x.lmr$lambdas[1]
      xT3[i] <- my.x.lmr$ratios[3]
      xT4[i] <- my.x.lmr$ratios[4]
      xT5[i] <- my.x.lmr$ratios[5]

      yL1[i] <- my.y.lmr$lambdas[1]
      yT2[i] <- my.y.lmr$lambdas[2]/my.y.lmr$lambdas[1]
      yT3[i] <- my.y.lmr$ratios[3]
      yT4[i] <- my.y.lmr$ratios[4]
      yT5[i] <- my.y.lmr$ratios[5]

      T2.12[i] <- mT2.12
      T2.21[i] <- mT2.21
      T3.12[i] <- mT3.12
      T3.21[i] <- mT3.21
      T4.12[i] <- mT4.12
      T4.21[i] <- mT4.21
      if(verbose) message("----------")
   }
   dev.off()
   if(verbose) message(" PDF X versus Y file=",xyfile," in ",getwd())


   if(verbose) message("Regional statistics and graphics commencing")
   ryL1 <- weighted.mean(yL1, weights=W)
   ryT2 <- weighted.mean(yT2, weights=W)
   ryT3 <- weighted.mean(yT3, weights=W)
   ryT4 <- weighted.mean(yT4, weights=W)
   ryT5 <- weighted.mean(yT5, weights=W)

   rxL1 <- weighted.mean(xL1, weights=W)
   rxT2 <- weighted.mean(xT2, weights=W)
   rxT3 <- weighted.mean(xT3, weights=W)
   rxT4 <- weighted.mean(xT4, weights=W)
   rxT5 <- weighted.mean(xT5, weights=W)

   rTAU <- weighted.mean(TAU, weights=W)
   rRHO <- weighted.mean(RHO, weights=W)

   rT2.12 <- weighted.mean(T2.12, weights=W)
   rT2.21 <- weighted.mean(T2.21, weights=W)
   rT3.12 <- weighted.mean(T3.12, weights=W)
   rT3.21 <- weighted.mean(T3.21, weights=W)
   rT4.12 <- weighted.mean(T4.12, weights=W)
   rT4.21 <- weighted.mean(T4.21, weights=W)
   if(verbose) message("    Regional weighted mean values of L-moments/L--comoments, done")

   z$sites <- each.site
   z$groups <- each.group
   z$record.lengths <- W
   z$kendall.tau <- TAU
   z$spearman.rho <- RHO

   z$lmoments <- list()

   z$lmoments$x <- list()
   z$lmoments$x$L1 <- xL1
   z$lmoments$x$T2 <- xT2
   z$lmoments$x$T3 <- xT3
   z$lmoments$x$T4 <- xT4
   z$lmoments$x$T5 <- xT5

   z$lmoments$y <- list()
   z$lmoments$y$L1 <- yL1
   z$lmoments$y$T2 <- yT2
   z$lmoments$y$T3 <- yT3
   z$lmoments$y$T4 <- yT4
   z$lmoments$y$T5 <- yT5

   z$lcomoments <- list()
   z$lcomoments$T2.12 <- T2.12
   z$lcomoments$T2.21 <- T2.21
   z$lcomoments$T3.12 <- T3.12
   z$lcomoments$T3.21 <- T3.21
   z$lcomoments$T4.12 <- T4.12
   z$lcomoments$T4.21 <- T4.21
   z$regional.means <- list()
   z$regional.means$x <- list()
   z$regional.means$x$L1 <- rxL1
   z$regional.means$x$T2 <- rxT2
   z$regional.means$x$T3 <- rxT3
   z$regional.means$x$T4 <- rxT4
   z$regional.means$x$T5 <- rxT5
   z$regional.means$y <- list()
   z$regional.means$y$L1 <- ryL1
   z$regional.means$y$T2 <- ryT2
   z$regional.means$y$T3 <- ryT3
   z$regional.means$y$T4 <- ryT4
   z$regional.means$y$T5 <- ryT5
   z$regional.means$kendall.tau <- rTAU
   z$regional.means$spearman.rho <- rRHO
   z$regional.means$T2.12 <- rT2.12
   z$regional.means$T2.21 <- rT2.21
   z$regional.means$T3.12 <- rT3.12
   z$regional.means$T3.21 <- rT3.21
   z$regional.means$T4.12 <- rT4.12
   z$regional.means$T4.21 <- rT4.21


   z$bivarsim <- NULL

   z$covariances <- list();

   xT2T3cov <- cov(cbind(xT2, xT3));
   xT3T4cov <- cov(cbind(xT3, xT4));
   xT2T4cov <- cov(cbind(xT2, xT4));

   yT2T3cov <- cov(cbind(yT2, yT3));
   yT3T4cov <- cov(cbind(yT3, yT4));
   yT2T4cov <- cov(cbind(yT2, yT4));
   z$covariances$x <- list();
   z$covariances$x$t2t3cov <- xT2T3cov;
   z$covariances$x$t3t4cov <- xT3T4cov;
   z$covariances$x$t2t4cov <- xT2T4cov;

   z$covariances$y <- list();
   z$covariances$y$t2t3cov <- yT2T3cov;
   z$covariances$y$t3t4cov <- yT3T4cov;
   z$covariances$y$t2t4cov <- yT2T4cov;
   if(verbose) message("    Covariances of the L-moments, done")


   T2T3cov.12 <- cov(cbind(T2.12, T3.12));
   T2T3cov.21 <- cov(cbind(T2.21, T3.21));
   T3T4cov.12 <- cov(cbind(T3.12, T4.12));
   T3T4cov.21 <- cov(cbind(T3.21, T4.21));
   T2T4cov.12 <- cov(cbind(T2.12, T4.12));
   T2T4cov.21 <- cov(cbind(T2.21, T4.21));


   z$covariances$t2t3cov.12 <- T2T3cov.12
   z$covariances$t2t3cov.21 <- T2T3cov.21
   z$covariances$t3t4cov.12 <- T3T4cov.12
   z$covariances$t3t4cov.21 <- T3T4cov.21
   z$covariances$t2t4cov.12 <- T2T4cov.12
   z$covariances$t2t4cov.21 <- T2T4cov.21

   if(verbose) message("    Covariances of the L-comoments, done")

   if(nsim > 0) {
     xT2T3bivn <- MASS::mvrnorm(nsim, mu = c(rxT2, rxT3),
                            Sigma = xT2T3cov)
     xT3T4bivn <- MASS::mvrnorm(nsim, mu = c(rxT3, rxT4),
                            Sigma = xT3T4cov)
     xT2T4bivn <- MASS::mvrnorm(nsim, mu = c(rxT2, rxT4),
                            Sigma = xT2T4cov)

     yT2T3bivn <- MASS::mvrnorm(nsim, mu = c(ryT2, ryT3),
                            Sigma = yT2T3cov)
     yT3T4bivn <- MASS::mvrnorm(nsim, mu = c(ryT3, ryT4),
                            Sigma = yT3T4cov)
     yT2T4bivn <- MASS::mvrnorm(nsim, mu = c(ryT2, ryT4),
                            Sigma = yT2T4cov)

     T2T3bivn.12 <- MASS::mvrnorm(nsim, mu = c(rT2.12, rT3.12),
                            Sigma = T2T3cov.12)
     T2T3bivn.21 <- MASS::mvrnorm(nsim, mu = c(rT2.21, rT3.21),
                            Sigma = T2T3cov.21)
     T3T4bivn.12 <- MASS::mvrnorm(nsim, mu = c(rT3.12, rT4.12),
                            Sigma = T3T4cov.12)
     T3T4bivn.21 <- MASS::mvrnorm(nsim, mu = c(rT3.21, rT4.21),
                            Sigma = T3T4cov.21)
     if(verbose) message("    Simulation of bivariate normals, done")
     z$bivarsim <- list()
     z$bivarsim$t2t3bivn.12 <- T2T3bivn.12
     z$bivarsim$t2t3bivn.21 <- T2T3bivn.21
     z$bivarsim$t3t4bivn.12 <- T3T4bivn.12
     z$bivarsim$t3t4bivn.21 <- T3T4bivn.21
     z$bivarsim$x <- list()
     z$bivarsim$x$t2t3 <- xT2T3bivn
     z$bivarsim$x$t3t4 <- xT3T4bivn
     z$bivarsim$x$t2t4 <- xT2T4bivn
     z$bivarsim$y <- list()
     z$bivarsim$y$t2t3 <- yT2T3bivn
     z$bivarsim$y$t3t4 <- yT3T4bivn
     z$bivarsim$y$t2t4 <- yT2T4bivn
   }

   if(verbose) message("*** lcoratdiagram(), done")

   return(z)
}

