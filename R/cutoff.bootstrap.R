cutoff.bootstrap<-function(xin,  M,  sim, dist, h.use, kfun, p1, p2, sig.lev )
{
  n.Boot <- boot(data= xin, statistic=S.n.Boot, R=M,  sim = sim,  dist=dist, h=h.use, kfun=kfun, p1=p1, p2=p2 )$t
  l.a.boot <- quantile(n.Boot, 1-sig.lev, type=7)
  l.a.boot
}
