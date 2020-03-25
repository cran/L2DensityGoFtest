\name{kde}
\title{Kernel Density Estimation}
\alias{kde}
\description{
Implements the (classical) kernel density estimator, see (2.2a) in  Silverman (1986).
}

\usage{
kde(xin, xout, h, kfun)
}

\arguments{
\item{xin}{A vector of data points.  Missing values not allowed.}
 \item{xout}{A vector of grid points at which the estimate will be calculated. }
 \item{h}{A scalar, the bandwidth to use in the estimate, e.g. \code{bw.nrd(xin)}}
 \item{kfun}{Kernel function to use. Supported kernels:  \code{\link{Epanechnikov}},  \code{\link{Biweight}},  \code{\link{Gaussian}},  \code{\link{Rectangular}},  \code{\link{Triangular}}. }
}

\details{
The classical kernel density estimator is given by
\deqn{\hat f(x;h) = n^{-1}\sum_{i=1}^n  K_h(x-X_{i})}
\eqn{h} is determined by  a bandwidth selector such as Silverman's default plug-in rule.
}


\value{
A vector with the density estimates at the designated points xout.
}

\references{
 Silverman (1986), Density Estimation for Statistics and Data Analysis, Chapman and Hall, London.

}

\author{
R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
}


\examples{
x<-seq(-5, 5,length=100)          #design points where the estimate will be calculated
plot(x, dnorm(x),  type="l", xlab = "x", ylab="density") #plot true density function
SampleSize <- 100
ti<- rnorm(SampleSize)            #draw a random sample from the actual distribution

huse<-bw.nrd(ti)
arg2<-kde(ti, x, huse, Epanechnikov) #Calculate the estimate
lines(x, arg2, lty=2)             #draw the result on the graphics device.
}
