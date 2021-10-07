\name{hopt.edgeworth}
\alias{hopt.edgeworth}


\title{Power-optimal bandwidth for density goodness-of-fit testing}

\description{
Implements the Edgeworth based power-optimal bandwidth for density goodness-of-fit testing.
}

\usage{
hopt.edgeworth(xin, dist, kfun, p1, p2, sig.lev)
}

\arguments{
\item{xin}{ A vector of data points - the available sample size.}
\item{dist}{The null distribution.}
\item{kfun}{The kernel to use in the density estimates used in the bandwidth expression.}
\item{p1}{Argument 1 (vector or object) for the null distribution.}
\item{p2}{Argument 2 (vector or object) for the null distribution.}
\item{sig.lev}{Significance level of the hypothesis test.}
}

\details{
Implements: the power-optimal bandwidth for the test statistic \code{\link{S.n}} given by
\deqn{ h = \left \{ \frac{\sqrt{2} K^{(3)}(0)}{3R(K)^{3/2}} \frac{\nu_2}{R(f)^{3/2}}\right \}^{-1/2} \left  \{ \frac{n \int \Delta_n^2 (x) f^2(x)\,dx}{\sigma^2 \{ 2 \nu_2 R(K)\}^{1/2}} \right \}^{-3/2}, }

which is similar in nature to the corresponding optimal bandwidth  value   obtained for the closely relatated regression setting in Gao and Gijbels (2008).

}

\value{
A scalar, the estimate the power-optimal bandwidth at the given significance level \eqn{\alpha}.
}

\references{
\href{https://amstat.tandfonline.com/doi/abs/10.1198/016214508000000968#.XnfOq-ozbCM}{Gao and Gijbels, Bandwidth selection in nonparametric kernel testing, pp. 1584-1594, JASA (2008)}
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
