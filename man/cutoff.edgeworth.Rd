\name{cutoff.edgeworth}
\alias{cutoff.edgeworth}


\title{Edgeworth based critical value for goodness-of-fit testing}

\description{
Implements the Edgeworth based critical value for goodness-of-fit testing of densities.
}

\usage{
cutoff.edgeworth(xin, dist, kfun, p1, p2, sig.lev)
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
Implements the Edgwworth based finite sample critical value, given by
\deqn{ l_\alpha = z_\alpha + d_0 \sqrt{h} + d_2(n \sqrt{h})^{-1} }
where \eqn{z_\alpha} is the \eqn{1-\alpha} quantile of the normal distribution and
\eqn{d_0 = d_1 - C_{ H_0}} and \eqn{d_j = (z_\alpha^2 - 1)c_j, j=1,2} with
\deqn{
c_1 = \frac{4K^{(3)}(0)\mu_2^3 \nu_3}{3\sigma^3}, \;  c_2 = \frac{\mu_3^2K^2(0)}{\sigma^3}, \; \mu_i =\int K^i(x)\,dx, i=1,\dots.
}
and
\deqn{
C_{H_0}   = 2\left (E f_0'( \theta_0) \right )^2 \Delta^{-1}, \; \nu_i   =  E \left \{f^{i}(x)\right \} = \int f^{i+1}(x)\,dx, i=1,\dots
}



This is similar in nature to the critical value estimate  obtained for the closely relatated regression setting in Gao and Gijbels (2008).
}


\value{
A scalar, the estimate the critical value at the given significance level
}

\references{
Bagkavos, Patil and Wood: Nonparametric goodness-of-fit testing for a continuous multivariate parametric model, (2021), under review.

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
