cutoff.edgeworth<-function(xin,   dist, kfun, p1, p2, sig.lev )
{
  n<-length(xin)
  h<-bw.nrd(xin) #pilot bandwidth for the kernel estimate

  DensityEst <- kde(xin,xin, h, Gaussian)
  Rhatfx <- mean(DensityEst) #estimate R(f)
  RhatfxUse<-Rhatfx^{3/2}

  RK<-  3/5 # for epanechnikov,  1/(2*sqrt(pi)) #for gaussian, 1/2 # for uniform,
  RKuse<-RK^{3/2}
  ConvK<- kernel.conv(kernel.conv(0, deriv.order = 0, kernel =  kfun)$kx, deriv.order = 0, kernel =  kfun)$kx #1/(2*sqrt(2*pi)) #fConvolution(kernel)
  ML.Dens<- NDistDens(xin, dist, p1, p2)
  mu.2 <- mean((DensityEst - ML.Dens)^2)
  nu.2 <- mean(DensityEst^2)
  SigmaSq <- 2* mu.2^2 *  nu.2 * RK

  banduse<-hopt.edgeworth(xin,   dist, kfun, p1, p2, sig.lev )
  RDeltaR<- mean( (ML.Dens - DensityEst)^2 * DensityEst )

  z.a<-qnorm(mean=mean(xin), sd=sd(xin), 1-sig.lev/2)
  n<-length(xin)
  sqrt.h<-sqrt(banduse)
  z.a.sq<-z.a^2 -1
  d0.1<- z.a.sq  *  ConvK * mean(DensityEst^3) *mu.2^3 / (3*SigmaSq^(3/2) )
  d0.2 <- RDeltaR /sqrt(2*nu.2 * RK)
  d0<- d0.1-d0.2
  d2<- z.a.sq *  mean((DensityEst - ML.Dens)^3)^2  * kernel.fun(0,0, kfun)$kx^2 / SigmaSq^(3/2)
  cut.off<- z.a + d0 * sqrt.h + d2/( n*sqrt.h)

  cut.off

}
