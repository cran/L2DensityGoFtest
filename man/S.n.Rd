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
\item{h}{The bandwidth to use, typically the output of \code{hopt.cutoff}.}
\item{drate}{power to use for the binning, the  default is 3/4.}
\item{dist}{The null distribution.}
\item{p1}{Argument 1 (vector or object) for the null distribution.}
\item{p2}{Argument 2 (vector or object) for the null distribution.}
}

\details{
Implements the test statistic used for testing the hypothesis
\deqn{H_0: f(x) = f_0(x, p1, p2) vs H_a: f(x) \neq f_0(x, p1, p2).}
This density goodness-of-fit test is based on a discretized approximation of the L2 distance. Assuming that \eqn{n} is the number of observations and \eqn{g = (max(xin)-min(xin))/n^{-drate}} is the number of bins in which the range of the data is split, the test statistic is:
\deqn{ S_n(h)  \frac{1}{g(g-1)h} {\sum\sum}_{i \neq j} K \{ (x_i-x_j)h^{-1}\} \{Y_i -f_0(x_i) \}\{Y_j -f_0(x_j) \} }
where \eqn{K} is the Epanechnikov kernel implemented in this package with the \code{\link{Epanechnikov}} function. The null model \eqn{f_0} is specified through the \code{dist} argument with parameters passed through the  \code{p1} and \code{p2} arguments.  The test is implemented in conjunction with \code{\link{hopt.cutoff}} function which provides the value of \eqn{h} needed for calculation of \eqn{ S_n(h)} and the critical value used to determine acceptance or rejection of the null hypothesis. See the example below for an application to a real world dataset.

The test statistic can be thought as a descritized version of the bias corrected test statistic in page 380 of Li and Racine (2007), using the nonparametric empirical density estimate in place of the classical kernel density estimate there.

}


\value{
A vector with the value of the test statistic as well as the bandwidth used for its calculation
}

\references{
\href{https://press.princeton.edu/books/hardcover/9780691121611/nonparametric-econometrics}{Li and Racine, Nonparametric Econometrics: Theory and Practice, Princeton University Press, (2007)}
}

\author{
R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
}

\seealso{ \code{\link{hopt.cutoff}}}

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
# exclude smallest value - only for uniform presentation of results (this step can be excluded)
xin = xin[xin!= min(xin)]

#inputs for the test statistic:
#kernel function to use in implementing the statistic and functional estimates for optimal h:
kfun<-"epanechnikov"
a.sig<-0.05 #define the significance level

Nulldist<-"normal" #null hypothesis is that the innovations are normaly distributed


p1<-mean(xin) # mean of the
p2<- sd(xin) #

cutoff<-hopt.cutoff(xin,   Nulldist, kfun, p1, p2, a.sig ) # power optimal bandwidth


TestStatistic<-S.n(xin, cutoff[1], 1/3, Nulldist, p1, p2)


cat("L2 test statistic value:", TestStatistic[1],  " critical value :", cutoff[2], "\n")
 # normality is rejected
}
