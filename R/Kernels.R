"Biweight"<-function(x, ...) { ifelse(abs(x) < 1, (15/16) * (  1 - x^2)^2, ifelse(abs(x) > 1, 0, 0)) }
"Epanechnikov"<-    function(x, ... ) { ifelse(abs(x)<1,  0.75 * (1-x^2) , ifelse(abs(x)>1, 0,0))  }
"Rectangular"<- function(x, ...){ ifelse(abs(x)<1, .5, ifelse(abs(x)>1,0,0))}
"Triangular"<-function(x, ...){ ifelse(abs(x)<1, 1-abs(x), ifelse(abs(x)>1,0,0)) }
"Gaussian"<- function(x, ...){ dnorm(x) }

