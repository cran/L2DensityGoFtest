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
Wand and Jones,  (1996), Kernel Smoothing, Chapman and Hall, London
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
