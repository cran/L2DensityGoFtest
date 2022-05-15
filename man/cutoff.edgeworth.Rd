\name{cutoff.edgeworth}
\alias{cutoff.edgeworth}


\title{Critical value based on Edgeworth expansion of the size function for the density goodness-of-fit test \eqn{\hat{S}_n(h)} of Bagkavos, Patil and Wood (2021)}

\description{
Implements the critical value for the density goodness-of-fit test \code{\link{S.n}}, approximating via an Edgeworth expansion the size function of the test statistic \code{\link{S.n}}.
}

\usage{
cutoff.edgeworth(xin, dist, kfun, p1, p2, sig.lev)
}

\arguments{
\item{xin}{ A vector of data points - the available sample.}
\item{dist}{The null distribution.}
\item{kfun}{The kernel to use in the density estimates used in the bandwidth expression.}
\item{p1}{Parameter 1 (vector or object) for the null distribution.}
\item{p2}{Parameter 2 (vector or object) for the null distribution.}
\item{sig.lev}{Significance level of the hypothesis test.}
}

\details{
Implements the critical value for the density goodness-of-fit test \code{\link{S.n}}, approximating via an Edgeworth expansion the size function of the test statistic \code{\link{S.n}}, given by
\deqn{ l_\alpha = z_\alpha + d_0 \sqrt{h} + d_2(n \sqrt{h})^{-1} }
where \eqn{z_\alpha} is the \eqn{1-\alpha} quantile of the normal distribution and
\eqn{d_0 = d_1 - C_{ H_0}} and \deqn{d_j = (z_\alpha^2 - 1)c_j, j=1,2} with
\deqn{
c_1 = \frac{4K^{(3)}(0)\mu_2^3 \nu_3}{3\sigma^3}, \;  c_2 = \frac{\mu_3^2K^2(0)}{\sigma^3}, \; \mu_i =\int K^i(x)\,dx, i=1,\dots.
}
and
\deqn{
C_{H_0}   = 2\left (E f_0'( \theta_0) \right )^2 \Delta^{-1}, \; \nu_i   =  E \left \{f^{i}(x)\right \} = \int f^{i+1}(x)\,dx, i=1,\dots
}
This critical value is the density function equivalent to the critical value estimate obtained in the closely relatated regression setting in Gao and Gijbels (2008) and is suitable for finite sample implementations of the test.
}


\value{
A scalar, the estimate of the critical value at the given significance level.
}

\references{
\href{https://amstat.tandfonline.com/doi/abs/10.1198/016214508000000968#.XnfOq-ozbCM}{Gao and Gijbels, Bandwidth selection in nonparametric kernel testing, pp. 1584-1594, JASA (2008)}
}


\author{
Dimitrios Bagkavos

R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
}


\seealso{
\code{\link{cutoff.asymptotic}, \link{cutoff.bootstrap} }
}


