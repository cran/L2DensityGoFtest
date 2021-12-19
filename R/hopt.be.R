hopt.be<-function(xin)
{
  n<-length(xin)
  h<-bw.SJ(xin) #pilot bandwidth for the kernel estimate
  DensityEst <- kde(xin,xin, h, Epanechnikov)
  Rhatfx <- mean(DensityEst) #estimate R(f)
  Rhatfxthird <- mean(DensityEst^3) #estimate R(f)
  banduse<- (1/sqrt(n))*sqrt((Rhatfx * 0.2571)/(Rhatfxthird * 0.345)) #based on the epanechnikov kernel
  banduse
}
