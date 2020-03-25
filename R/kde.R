"kde"<-function(xin, xout, h, kfun)
{
  n<- length(xin)
  arg1<-(sapply(xout, "-", xin))/h  
  arg2<-kfun(arg1)
  arg3<-colSums(arg2)
  arg3/(n*h)
}
