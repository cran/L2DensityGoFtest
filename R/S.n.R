S.n<-function(xin, h,  dist, p1, p2)
  #xin: the random sample
  # h: the bandwidth to use in the test statistic - also conider the very first MSE optimal rule
  # dist: distribution under the null
  # p1: 1st parameter of the distribution dist. In case of Normal Mixtures this is the object MW.nm7 etc
  # p2: 2nd parameter of the distribution dist.
{
  MinX<-min(xin)
  MaxX<-max(xin)
  N<-length(xin)
  Delta<-1/N^{3/4} #n is the number of bins - set by A1.
  n<- ceiling( (MaxX - MinX)/Delta )

  PartInt<-c(MinX, MinX+(1:n)*Delta) #Partion the interval (x[1], x[n]) into subintervals
  BinCenters<-( PartInt[1:n] + PartInt[2:(n+1)])/2 # calculate the Bincenters:
  yi <- sapply(1:n, function(i, xin, PartInt) length(which(xin >= PartInt[i] & xin < PartInt[i+1])), xin, PartInt)/(N*Delta) #scaled bin counts as p.d.f. estimates
  arg1<-(sapply(BinCenters, "-", BinCenters))  # cross product of all BinCenters in order to calculate K{(x_i - x_j)/h}
  Null.Dens.Est<- NDistDens(BinCenters, dist,   p1, p2)
  Dens.Diffs<- yi-Null.Dens.Est  # Yi - f(x_i)
  All.Dens.Diffs<- Dens.Diffs %*% t(Dens.Diffs)  # cross product of all Yi - f(x_i) diffs
  arg2<-Epanechnikov(arg1/h)
  test.stat.tmp<-arg2 * All.Dens.Diffs

  test.statistic<- sum(test.stat.tmp)
  c( test.statistic, Delta)
}
