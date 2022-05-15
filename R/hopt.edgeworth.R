hopt.edgeworth<-function(xin,   dist, kfun, p1, p2, sig.lev )
{
  n<-length(xin)
  h<-bw.nrd(xin) #pilot bandwidth for the kernel estimate

  DensityEst <- kde(xin,xin, h, Gaussian)
  Rhatfx <- mean(DensityEst) #estimate R(f)
  RhatfxUse<-Rhatfx^{3/2}

  RK<-  3/5 # for epanechnikov,  1/(2*sqrt(pi)) #for gaussian, 1/2 # for uniform,
  RKuse<-RK^{3/2}
  ConvK<- switch(kfun,
                 Gaussian = 0.2765382,
                 Epanechnikov = 0.409542,
                 Triangular = 0.3703704,
                 Rectangular = 0.3703704,
                 Biweight = 0.3276294,
                 Epanechnikov2 = 0.409542
                 ) #kernel.conv(kernel.conv(0, deriv.order = 0, kernel =  kfun)$kx, deriv.order = 0, kernel =  kfun)$kx
  ML.Dens<- NDistDens(xin, dist, p1, p2)
  mu.2 <- mean((DensityEst - ML.Dens)^2)
  nu.2 <- mean(DensityEst^2)
  SigmaSq <- 2* mu.2^2 *  nu.2 * RK # 16/45 * RK * Rhatfx #estimate test statistic variance

  term1 <- ( (sqrt(2)*ConvK*nu.2)/(3*RKuse*RhatfxUse) )^{-1/2} #
  RDeltaR<- mean( (ML.Dens - DensityEst)^2 * DensityEst )
  term2<- ( 1/(    sqrt(2*   nu.2 * RK)) )^{-3/2}
  banduse<-   term1 * (n*RDeltaR)^{-3/2} *term2

  banduse

}
