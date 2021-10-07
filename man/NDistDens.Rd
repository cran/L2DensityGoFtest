\name{NDistDens}
\alias{NDistDens}


\title{Select null distribution}

\description{
Implements the selection of null distribution; to be used within the implementation of the test statistic S.n
}

\usage{
NDistDens(x, dist, p1, p2)
}

\arguments{
\item{x}{ A vector of data points - the available sample size.}
\item{dist}{The null distribution.}
\item{p1}{Argument 1 (vector or object) for the null distribution.}
\item{p2}{Argument 2 (vector or object) for the null distribution.}
}

\details{
Implements the null distribution evaluation at designated points, given the parameters p1 and p2.
}

\value{
A vector containing the density values of the designated distribution
}

\references{
Bagkavos, Patil and Wood: Nonparametric goodness-of-fit testing for a continuous multivariate parametric model, (2021), under review.
}

\author{
Dimitrios Bagkavos

R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
}

%\note{
%%  ~~further notes~~
%}

%\seealso{
%\code{\link{alfa.profile}, \link{alfa}, \link{alfainv}, \link{alpha.mle}
%}
%}

%\examples{
%
%}
