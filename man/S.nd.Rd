\name{S.nd}
\title{Density goodness-of-fit test statistic based on discretized L2 distance}
\alias{S.nd}

\description{
Implements the multivariate (d >=2) density goodness of fit test statistic \eqn{\hat{S}_n(h)} of Bagkavos, Patil and Wood (2021), based on aggregation of local discrepancies between the fitted parametric density and a nonparametric empirical density estimator.
}

\usage{
S.nd(xin, h,  dist, p1, p2)
}

\arguments{
\item{xin}{A matrix (n x d) of data points - the available sample with n rows and d columns, each column corresponds to a different coordinate axis.}
\item{h}{The bandwidth vector to use, typically the output of \code{hopt.be} in each coordinate direction.}
\item{dist}{The null distribution.}
\item{p1}{Parameter 1 (vector or object) for the null distribution.}
\item{p2}{Parameter 2 (vector or object) for the null distribution.}
}

\details{
Implements the test statistic used for testing the hypothesis
\deqn{H_0: f(x) = f_0(x, p1, p2) \;\; vs \;\; H_a: f(x) \neq f_0(x, p1, p2).}
This density goodness-of-fit test is based on a discretized approximation of the L2 distance. Assuming that \eqn{n} is the number of observations and \eqn{g = (max(xin)-min(xin))/n^{-drate}} is the number of bins in which the range of the data is split, the test statistic is:
\deqn{ S_n(h) = n \Delta^2  {\sum\sum}_{i \neq j} K \{ (X_{i1}-X_{j1})h_1^{-1}, \dots,  (X_{id}-X_{jd})h_d^{-1} \} \{Y_i -f_0(X_i) \}\{Y_j -f_0(X_j) \} }
where \eqn{K} is the Epanechnikov kernel implemented in this package with the \code{\link{Epanechnikov}} function. The null model \eqn{f_0} is specified through the \code{dist} argument with parameters passed through the  \code{p1} and \code{p2} arguments.  The test is implemented either with bandwidth \code{\link{hopt.edgeworth}} or with  bandwidth \code{\link{hopt.be}} which provide the value of \eqn{h} needed for calculation of \eqn{ S_n(h)} and the critical value used to determine acceptance or rejection of the null hypothesis.
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

\seealso{ \code{\link{S.n}}}

\examples{
library(mvtnorm)
sigma <- matrix(c(4,2,2,3), ncol=2)

x <- rmvnorm(n=100, mean=c(1,2), sigma=sigma)
h.be1 <- hopt.be(x[,1])
h.be2 <- hopt.be(x[,2])
h<-c(h.be1, h.be2)
Nulldist<-"normal"

S.nd(x, h,  Nulldist, c(1,2), sigma)
}
