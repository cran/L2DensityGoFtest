\name{S.n}
\title{Density goodness-of-fit test statistic based on discretized L2 distance}
\alias{S.n}

\description{
Implements the density goodness of fit test statistic \eqn{\hat{S}_n(h)} of Bagkavos, Patil and Wood (2021), based on aggregation of local discrepancies between the fitted parametric density and a nonparametric empirical density estimator.
}

\usage{
S.n(xin, h,  dist, p1, p2)
}

\arguments{
\item{xin}{A vector of data points - the available sample size.}
\item{h}{The bandwidth to use, typically the output of \code{hopt.edgeworth}.}
\item{dist}{The null distribution.}
\item{p1}{Parameter 1 (vector or object) for the null distribution.}
\item{p2}{Parameter 2 (vector or object) for the null distribution.}
}

\details{
Implements the test statistic used for testing the hypothesis
\deqn{H_0: f(x) = f_0(x, p1, p2) \;\; vs \;\; H_a: f(x) \neq f_0(x, p1, p2).}
This density goodness-of-fit test is based on a discretized approximation of the L2 distance. Assuming that \eqn{n} is the number of observations and \eqn{g = (max(xin)-min(xin))/n^{-drate}} is the number of bins in which the range of the data is split, the test statistic is:
\deqn{ S_n(h) = n \Delta^2 h^{-1/2} {\sum\sum}_{i \neq j} K \{ (X_i-X_j)h^{-1}\} \{Y_i -f_0(X_i) \}\{Y_j -f_0(X_j) \} }
where \eqn{K} is the Epanechnikov kernel implemented in this package with the \code{\link{Epanechnikov}} function. The null model \eqn{f_0} is specified through the \code{dist} argument with parameters passed through the  \code{p1} and \code{p2} arguments.  The test is implemented either with bandwidth \code{\link{hopt.edgeworth}} or with  bandwidth \code{\link{hopt.be}} which provide the value of \eqn{h} needed for calculation of \eqn{ S_n(h)} and the critical value used to determine acceptance or rejection of the null hypothesis. See the example below for an application to a real world dataset.

}


\value{
A vector with the value of the test statistic as well as the Delta value used for its calculation
}

\references{
Bagkavos, Patil and Wood: Nonparametric goodness-of-fit testing for a continuous multivariate parametric model, (2021), under review.
}

\author{
R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
}

\seealso{ \code{\link{S.n.Boot}}}

\examples{
library(fGarch)
library(kedd)
library(boot)
 \dontrun{data(EuStockMarkets)
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
# Edgeworth cutoff point:
cutoff<-cutoff.edgeworth(xin,   Nulldist, kfun, p1, p2, a.sig )
# Bootstrap cutoff point:
cutoff.boot<-cutoff.bootstrap(xin, 100,  "permutation", Nulldist, h.be, kfun, p1, p2, a.sig)
# Asympt. Norm. cutoff point:
cutoff.asympt<-cutoff.asymptotic( Nulldist,   p1, p2, a.sig )

TestStatistic<-S.n(xin, h, Nulldist, p1, p2)
TestStatistic.be<-S.n(xin, h.be, Nulldist, p1, p2)

cat("L2 test statistic value with power opt. band:", TestStatistic[1],
"\nL2 test statistic value Barry-Essen bandwidth:", TestStatistic.be[1],
"\ncritical value asymptotic:", round(cutoff.asympt,3), "critical value bootstrap:",
round(cutoff.boot,3),  "critical value Edgeworth:", round(cutoff,3), "\n")
#L2 test statistic value Edgeworth: 7.257444
#L2 test statistic value Berry-Esseen bandwidth: 10.97069
# critical value Asymptotically Norm.:  1.801847
# critical value Edgeworth: 2.140446
# critical value bootstrap: 6.040048
# L2 test statistic >  critical value on all occasions, hence normality is rejected}
}
