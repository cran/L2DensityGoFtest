\name{Kernels}
\alias{Biweight}
\alias{Epanechnikov}
\alias{Triangular}
\alias{Gaussian}
\alias{Rectangular}

\title{Kernel functions}

\description{
Implements various kernel functions, including boundary, integrated and discrete kernels for use in the definition of the nonparametric estimates
}

\usage{
Biweight(x, ...)
Epanechnikov(x, ...)
Triangular(x, ...)
Gaussian(x, ...)
Rectangular(x, ...)
}

\arguments{
\item{x}{ A vector of data points where the kernel will be evaluated.}
\item{...}{Further arguments.}
}

\details{
Implements the Biweight, Triangular, Guassian, Rectangular and Epanechnikov kernels.
}

\value{
The value of the kernel at \eqn{x}
}

\references{
\enumerate{
\item \href{https://ieeexplore.ieee.org/document/4385743}{Bagkavos and Patil, Local Polynomial Fitting in Failure Rate Estimation, IEEE Transactions on Reliability, 57, (2008)},
\item \href{https://link.springer.com/article/10.1007/s10463-010-0277-6}{Bagkavos (2011),  Annals of the Institute of Statistical Mathematics, 63(5), 1019-1046},
}
}

%\author{
%Dimitrios Bagkavos
%
%R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
%}

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
