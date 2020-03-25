NDistDens<-function(x, dist, p1, p2)
{
  switch(dist,
         normixt = dnorMix(x,  p1),
         exponential = dexp(x, p1),
         lnorm = dlnorm(x,  p1, p2),
         weibull = dweibull(x, p1, p2),
         normal = dnorm(x, p1, p2)
  )
}
