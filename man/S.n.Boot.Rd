\name{S.n.Boot}
\title{Goodness-of-Fit test statistic based on discretized L2 distance}
\alias{S.n.Boot}

\description{
Implements the bootstraped version of the density goodness-of-fit test \eqn{\hat{S}_n(h)} defined in (6) Bagkavos, Patil and Wood (2021).
}

\usage{
S.n.Boot(xin1, indices, h,  dist, kfun, p1, p2)
}

\arguments{
\item{xin1}{A vector of data points to perfrom bootstrap on.}
\item{indices}{indices to use for the bootstrap process.}
\item{h}{The bandwidth to use, typically the output of \code{hopt.be}.}
\item{dist}{The null distribution.}
\item{kfun}{The kernel to use in the density estimates used in the bandwidth expression.}
\item{p1}{Argument 1 (vector or object) for the null distribution.}
\item{p2}{Argument 2 (vector or object) for the null distribution.}
}

\details{
Implements the bootstrap version of the test statistic \code{S.n} for use in the  \code{\link{cutoff.bootstrap}} function. This function is typically not to be called directly by the user; it is rather meant to be called indirectly through the \code{cutoff.bootstrap} function.
}


\value{
A vector of values of the test statistic.
}

\references{
Bagkavos, Patil and Wood: Nonparametric goodness-of-fit testing for a continuous multivariate parametric model, (2021), under review.
}

\author{
R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
}

\seealso{ \code{\link{S.n}}}


