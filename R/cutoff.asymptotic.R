cutoff.asymptotic<-function( dist,   p1, p2, sig.lev )
{
  RK<-  3/5
  asymptotic.cut.off<- 2*RK* (integrate(NDistDens, lower=-Inf, upper=Inf, dist=dist, p1=p1, p2=p2)$value)^2 *qnorm(1-sig.lev)
  asymptotic.cut.off
}
