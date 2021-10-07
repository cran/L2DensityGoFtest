\name{S.n}
\title{Goodness-of-Fit test statistic based on discretized L2 distance}
\alias{S.n}

\description{
Implements a density goodness of fit test based on a discretized approximation of the L2 distance.
}

\usage{
S.n(xin, h, drate, dist, p1, p2)
}

\arguments{
\item{xin}{A vector of data points - the available sample size.}
\item{h}{The bandwidth to use, typically the output of \code{hopt.edgeworth}.}
\item{drate}{power to use for the binning, the  default is 3/4.}
\item{dist}{The null distribution.}
\item{p1}{Argument 1 (vector or object) for the null distribution.}
\item{p2}{Argument 2 (vector or object) for the null distribution.}
}

\details{
Implements the test statistic used for testing the hypothesis
\deqn{H_0: f(x) = f_0(x, p1, p2) \;\; vs \;\; H_a: f(x) \neq f_0(x, p1, p2).}
This density goodness-of-fit test is based on a discretized approximation of the L2 distance. Assuming that \eqn{n} is the number of observations and \eqn{g = (max(xin)-min(xin))/n^{-drate}} is the number of bins in which the range of the data is split, the test statistic is:
\deqn{ S_n(h) = n \Delta^2 h^{-1/2} {\sum\sum}_{i \neq j} K \{ (X_i-X_j)h^{-1}\} \{Y_i -f_0(X_i) \}\{Y_j -f_0(X_j) \} }
where \eqn{K} is the Epanechnikov kernel implemented in this package with the \code{\link{Epanechnikov}} function. The null model \eqn{f_0} is specified through the \code{dist} argument with parameters passed through the  \code{p1} and \code{p2} arguments.  The test is implemented either with bandwidth \code{\link{hopt.edgeworth}} or with  bandwidth \code{\link{hopt.be}} which provide the value of \eqn{h} needed for calculation of \eqn{ S_n(h)} and the critical value used to determine acceptance or rejection of the null hypothesis. See the example below for an application to a real world dataset.

}


\value{
A vector with the value of the test statistic as well as the bandwidth used for its calculation
}

\references{
Bagkavos, Patil and Wood: Nonparametric goodness-of-fit testing for a continuous multivariate parametric model, (2021), under review.
}

\author{
R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
}

\seealso{ \code{\link{hopt.edgeworth}}}

\examples{
library(fGarch)
library(kedd)
data(EuStockMarkets)
DAX <- as.ts(EuStockMarkets[,"DAX"])
dax <-  diff(log(DAX))#[,"DAX"]

# Fit a GARCH(1,1) model to dax returns:
lll<-garchFit(~ garch(1,1), data = as.ts(dax), trace = FALSE, cond.dist ="std")
# define the model innovations, to be used as input to the test statistic
xin<-lll@residuals /lll@sigma.t
# exclude smallest value - only for uniform presentation of results
#(this step can be excluded):
xin = xin[xin!= min(xin)]

#inputs for the test statistic:
#kernel function to use in implementing the statistic
#and functional estimates for optimal h:
kfun<-"epanechnikov"
a.sig<-0.05 #define the significance level
#null hypothesis is that the innovations are normaly distributed:
Nulldist<-"normal"

p1<-mean(xin)
p2<- sd(xin)
#Power optimal bandwidth:
h<-hopt.edgeworth(xin,   Nulldist, kfun, p1, p2, a.sig )
h.be <- hopt.be(xin)
# Edgeworth based cutoff point:
cutoff<-cutoff.edgeworth(xin,   Nulldist, kfun, p1, p2, a.sig )

TestStatistic<-S.n(xin, h, 1/3, Nulldist, p1, p2)
TestStatistic.be<-S.n(xin, h.be, 1/3, Nulldist, p1, p2)

cat("L2 test statistic value Edgeworth:", TestStatistic[1],
"\nL2 test statistic value Barry-Essen bandwidth:", TestStatistic.be[1],
"\ncritical value :", cutoff, "\n")
#L2 test statistic value Edgeworth: 21.49267
#L2 test statistic value Barry-Essen bandwidth: 52.66799
# critical value : 2.140446
# L2 test statistic >  critical value, hence normality is rejected
}
