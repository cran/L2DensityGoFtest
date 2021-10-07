\name{cutoff.asymptotic}
\alias{cutoff.asymptotic}


\title{Asymptotic critical value for goodness-of-fit testing}

\description{
Implements an asymptotic based critical value for goodness-of-fit testing of densities.
}

\usage{
cutoff.asymptotic(dist,  p1, p2, sig.lev)
}

\arguments{
\item{dist}{The null distribution.}
\item{p1}{Argument 1 (vector or object) for the null distribution.}
\item{p2}{Argument 2 (vector or object) for the null distribution.}
\item{sig.lev}{Significance level of the hypothesis test.}
}

\details{
Implements the asymptotic critical value in Remark 1, Bagkavos, Patil and Wood (2021), equal to  \eqn{ z_\alpha \sigma_{0, \theta_0} } where \eqn{z_\alpha} is the \eqn{1-\alpha} quantile of the normal distribution and
\deqn{ \sigma_{0, \theta_0}^2 = 2 \left (\int K^2(u)\,du \right ) \left (\int f^2_0(x; \theta_0)\,dx \right ).  }

}


\value{
A scalar, the estimate the critical value at the given significance level
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
