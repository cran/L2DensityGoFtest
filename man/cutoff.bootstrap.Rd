\name{cutoff.bootstrap}
\alias{cutoff.bootstrap}


\title{Bootstrap critical value for the goodness-of-fit test statistic \eqn{\hat{S}_n(h)} of Bagkavos, Patil and Wood (2021)}

\description{
Implements a bootstrap critical value for testing the goodness-of-fit of a parametrically estimated density with the test statistic \code{S.n}.
}

\usage{
cutoff.bootstrap(xin, M,  sim, dist, h.use, kfun, p1, p2, sig.lev)
}


\arguments{
\item{xin}{ A vector of data points - the available sample.}
\item{M}{Number of bootstrap replications.}
\item{sim}{A character string indicating the type of simulation required: "ordinary" (the default), "parametric", "balanced", "permutation", or "antithetic".}
\item{dist}{The null distribution.}
\item{h.use}{The test statistic bandwidth, best implemented with \code{\link{hopt.be}}.}
\item{kfun}{The kernel to use in the density estimates used in the bandwidth expression.}
\item{p1}{Parameter 1 (vector or object) for the null distribution.}
\item{p2}{Parameter 2 (vector or object) for the null distribution.}
\item{sig.lev}{Significance level of the hypothesis test.}
}

\details{
Implements the bootstrap based finite sample critical value defined in Section 2.6, Bagkavos, Patil and Wood (2021), and calculated as follows:

1.  Resample the observations \eqn{\mathcal{X}=\{X_1, \dots, X_n\}} to obtain \eqn{M} bootstrap samples, denoted by \eqn{\mathcal{X}_m^\ast=\{ X_{1m}^\ast, \dots, X_{nm}^\ast\}},  where for each \eqn{m=1,\ldots , M}, \eqn{\mathcal{X}_m^\ast} is sampled randomly, with replacement, from \eqn{\mathcal{X}}.  Write \eqn{\hat{\theta}=\theta(\mathcal{X})} for the estimator of \eqn{\theta} based on the original sample \eqn{\mathcal{X}} and, for each \eqn{m}, define the bootstrap estimator of \eqn{\theta} by \eqn{\hat{\theta}_m^\ast = \theta(\mathcal{X}_m^\ast)}, where \eqn{\theta(\cdot)} is the relevant functional for the parameter \eqn{\theta}.

2. For \eqn{m=1, \ldots , M}, use \eqn{\mathcal{X}_m^\ast =\{X_{1m}^\ast, \dots, X_{nm}^\ast\}} and \eqn{\hat \theta_m^\ast}  from the previous step to calculate   \eqn{n \Delta^{2d} h^{-d/2} \hat S_{n,m}^\ast(h\rho)},\eqn{m=1, \dots, M}.

3. Calculate   \eqn{\ell_\alpha^\ast}  as the \eqn{1-\alpha} empirical quantile of the values \eqn{n \Delta^{2d} h^{-d/2} \hat S_{n,m}^\ast(h\rho)}, \eqn{m=1, \dots, M}. Then \eqn{\ell_\alpha^\ast}   approximately   satisfies \eqn{P^\ast [ n  \Delta^{2d} h^{-d/2}\hat S_{n,m}^\ast(h\rho)> \ell_\alpha^\ast ]=1-\alpha}, where \eqn{P^\ast} indicates the bootstrap probability measure conditional on \eqn{\mathcal{X}}.

}





\value{
A scalar, the estimate of the bootstrap critical value at the given significance level.
}

\references{
Bagkavos, Patil and Wood: Nonparametric goodness-of-fit testing for a continuous multivariate parametric model, (2021), under review.

\href{https://amstat.tandfonline.com/doi/abs/10.1198/016214508000000968#.XnfOq-ozbCM}{Gao and Gijbels, Bandwidth selection in nonparametric kernel testing, pp. 1584-1594, JASA (2008)}
}


\author{
Dimitrios Bagkavos

R implementation and documentation: Dimitrios Bagkavos <dimitrios.bagkavos@gmail.com>
}



\seealso{
\code{\link{cutoff.asymptotic}, \link{cutoff.edgeworth} }
}

\examples{
library(nor1mix)
library(boot)
SampleSize<-80
M<-1000
dist<- "normixt"
kfun<- Epanechnikov
p1 <-MW.nm2
p2 <-1
sig.lev <- 0.05

sim<-"ordinary"

#Run the following to compare the asymptotic and bootstrap cut-off points on 4 occasions:
for(i in 15:18)
  {
    set.seed(i)
    xin<-rnorMix(SampleSize, p1)
    h.use <- hopt.be(xin)
    l.a.a<-cutoff.asymptotic( dist,   p1, p2, sig.lev )
    l.a.b<- cutoff.bootstrap(xin,  M,  sim, dist, h.use,  kfun, p1, p2, sig.lev)
    #print the result of each iteration:
    cat("Asympt. cut.off= ", l.a.a, "Boot. cut.off= ", l.a.b,  "\n")
   }
}
