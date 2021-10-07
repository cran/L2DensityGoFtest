\name{hopt.be}
\alias{hopt.be}


\title{Power-optimal bandwidth and critical value}

\description{
Implements the power-optimal bandwidth for density goodness-of-fit testing based on a Berry Essen bound.
}

\usage{
hopt.be(xin)
}

\arguments{
\item{xin}{ A vector of data points - the available sample size.}
}

\details{
Implements the power-optimal bandwidth for the test statistic \code{\link{S.n}} given by
\deqn{ h =   n^{-1/2} \sqrt{\frac{\hat \nu_p R_4(K)}{\rho_\ast^2 \hat \nu_4 I_0(K)} },   }
where
\deqn{  \hat \nu_p = n^{-1} \sum_{j=1}^n \hat f(X_j; \hat h_a),  }

and \eqn{\hat h_a} is the density optimal bandwidth calculated by a reference to a prametric distribution, \eqn{\rho_\star=1} and
\deqn{   R_4(K)=\int K^4(x)\,dx.}
}
\value{
The estimate of the power-optimal bandwidth
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
