\name{hopt.be}
\alias{hopt.be}


\title{Power-optimal bandwidth for the test statistic \eqn{\hat{S}_n(h)}}

\description{
Implements an optimal, with respect to Berry-Esseen bound, bandwidth for the density goodness-of-fit test \eqn{\hat{S}_n(h)} of Bagkavos, Patil and Wood (2021).
}

\usage{
hopt.be(xin)
}

\arguments{
\item{xin}{ A vector of data points - the available sample.}
}

\details{
Implements the  Berry-Esseen bound optimal bandwidth  defined in (18), Bagkavos, Patil and Wood (2022), given by
\deqn{ h =   n^{-1/2} \sqrt{\frac{\hat \nu_p R_4(K)}{\rho_\ast^2 \hat \nu_4 I_0(K)} },   }
where
\deqn{  \hat \nu_p = n^{-1} \sum_{j=1}^n \hat f(X_j; \hat h_a),  }

and \eqn{\hat h_a} is the density optimal bandwidth calculated by a reference to a prametric distribution, \eqn{\rho_\star=1} and
\deqn{   R_4(K)=\int K^4(x)\,dx.}
}
\value{
The estimate of the Berry-Esseen optimal bandwidth.
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

\seealso{
\code{\link{hopt.edgeworth}  }
}

%\examples{
%
%}
