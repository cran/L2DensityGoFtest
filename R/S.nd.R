S.nd<-function(xin, h,  dist, p1, p2)
{
  d<-ncol(xin)
  res<- 1
  n.use<-vector("numeric", length=d)
  for(i in 1:d)
  {
    xin.use<-xin[,i]
    MinX<-min(xin.use)
    MaxX<-max(xin.use)
    N<-length(xin.use)
    Delta<-1/N^{3/4} #n is the number of bins - set by A1.
    n<- ceiling( (MaxX - MinX)/Delta )
    n.use[i]<-n
  }
  n<-max(n.use)
  for(i in 1:d)
  {
    xin.use<-xin[,i]
    h.use<-h[i]
    MinX<-min(xin.use)
    MaxX<-max(xin.use)
    N<-length(xin.use)
    Delta<-1/N^{3/4} #n is the number of bins - set by A1.

    PartInt<-c(MinX, MinX+(1:n)*Delta) #Partion the interval (x[1], x[n]) into subintervals
    BinCenters<-( PartInt[1:n] + PartInt[2:(n+1)])/2 # calculate the Bincenters:
    yi <- sapply(1:n, function(i, xin.use, PartInt) length(which(xin.use >= PartInt[i] & xin.use < PartInt[i+1])), xin.use, PartInt)/(N*Delta) #scaled bin counts as p.d.f. estimates
    cat(length(yi),"\n")
    arg1<-(sapply(BinCenters, "-", BinCenters))  # cross product of all BinCenters in order to calculate K{(x_i - x_j)/h}
    Null.Dens.Est<- NDistDens(BinCenters, dist,   p1[i], p2[i,i])
    Dens.Diffs<- yi-Null.Dens.Est  # Yi - f(x_i)

     All.Dens.Diffs<- Dens.Diffs %*% t(Dens.Diffs)  # cross product of all Yi - f(x_i) diffs
     arg2<-Epanechnikov(arg1/h.use)
     test.stat.tmp<-arg2   * All.Dens.Diffs
     res <- res * test.stat.tmp  * h.use^{-1/2}

  }
  result_final <- sum(  res)
  test.statistic<-  result_final * N * Delta^2
   c( test.statistic, Delta)
}

