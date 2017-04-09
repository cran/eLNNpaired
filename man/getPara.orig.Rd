\name{getPara.orig}
\alias{getPara.orig}
\title{
Translate Re-Parameterized Parameters to Original Scale
}
\description{
Translate re-parameterized parameters to original scale.
}
\usage{
getPara.orig(
  delta1, 
  xi1, 
  lambda1, 
  nu1, 
  delta2, 
  xi2, 
  lambda2, 
  nu2, 
  lambda3, 
  nu3, 
  c1 = qnorm(0.95), 
  c2 = qnorm(0.05))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{delta1}{
log of the mean of the mean expression levels for gene probes in cluster 1 (over-expressed probes).
}
  \item{xi1}{
the value of the inverse function of the cumulative distribution
function of the standard normal distribution at
the point that is equal to the scalar in the variance of the mean expression levles for gene probes
in cluster 1 (over-expressed probes).
}
  \item{lambda1}{
a parameter related to \eqn{\alpha_1}{alpha_1}, which is
the shape parameter of the distribution of the variance of gene expression
levels for gene probes in cluster 1 (over-expressed probes).
}
  \item{nu1}{
log of the rate parameter of the distribution of the variance of gene expression
levels for gene probes in cluster 1 (over-expressed probes).
}
  \item{delta2}{
log of the negative mean of the mean expression levels for gene probes in cluster 2 (under-expressed probes).

}
  \item{xi2}{
the value of the inverse function of the cumulative distribution
function of the standard normal distribution at
the point that is equal to the scalar in the variance of the mean expression levles for gene probes
in cluster 2 (under-expressed probes).
}
  \item{lambda2}{
a parameter related to \eqn{\alpha_2}{alpha_2}, which is
the shape parameter of the distribution of the variance of gene expression
levels for gene probes in cluster 2 (under-expressed probes).
}
  \item{nu2}{
log of the rate parameter of the distribution of the variance of gene expression
levels for gene probes in cluster 2 (under-expressed probes).

}
  \item{lambda3}{
a parameter related to \eqn{\alpha_3}{alpha_3}, which is
the shape parameter of the distribution of the variance of gene expression
levels for gene probes in cluster 3 (non-differentially-expressed probes).
}
  \item{nu3}{
log of the rate parameter of the distribution of the variance of gene expression
levels for gene probes in cluster 3 (non-differentially-expressed probes).
}
  \item{c1}{
the lower bound for \eqn{\mu_g/\sqrt{\tau_g^{-1}}}{mu_g/sqrt(tau_g^{-1})}
for cluster 1 (over-expressed probes).
By default \eqn{c_1=\Phi^{-1}(0.95)}{c_1=Phi^{-1}(0.95)}.
}
  \item{c2}{
the upper bound for \eqn{\mu_g/\sqrt{\tau_g^{-1}}}{mu_g/sqrt(tau_g^{-1})}
for cluster 2 (under-expressed probes).
By default \eqn{c_2=\Phi^{-1}(0.05)}{c_2=Phi^{-1}(0.05)}.
}
}
\details{
We assume the following the Bayesian hierarchical models 
for the 3 clusters of gene probes.

For cluster 1 (over-expressed gene probes):
\deqn{d_{gl}|\left(\mu_g, \tau_g\right)\sim
N\left(\mu_g, \tau_g^{-1}\right),\\
\mu_g | \tau_g \sim N\left(\mu_1, k_1 \tau_g^{-1}\right),\\
\tau_g\sim \Gamma\left(\alpha_1, \beta_1\right).
}{
d_{gl} | (mu_g, tau_g) ~
N(mu_g, tau_g^{-1}),\\
mu_g | tau_g ~ N(mu_1, k_1 tau_g^{-1}),\\
tau_g ~ Gamma(alpha_1, beta_1).
}

For cluster 2 (under-expressed gene probes):
\deqn{d_{gl}|\left(\mu_g, \tau_g\right)\sim
N\left(\mu_g, \tau_g^{-1}\right),\\
\mu_g | \tau_g \sim N\left(\mu_2, k_2 \tau_g^{-1}\right),\\
\tau_g\sim \Gamma\left(\alpha_2, \beta_2\right).
}{
d_{gl} | (mu_g, tau_g) ~
N(mu_g, tau_g^{-1}),\\
mu_g | tau_g ~ N(mu_2, k_2 tau_g^{-1}),\\
tau_g ~ Gamma(alpha_2, beta_2).
}

For cluster 3 (non-differentially-expressed gene probes):
\deqn{d_{gl}|\left(\tau_g\right)\sim
N\left(0, \tau_g^{-1}\right),\\
\tau_g\sim \Gamma\left(\alpha_3, \beta_3\right).
}{
d_{gl} | (mu_g, tau_g) ~
N(0m, tau_g^{-1}),\\
tau_g ~ Gamma(alpha_3, beta_3).
}

For cluster 1, we add one constraint
\deqn{
\alpha_1>1+\beta_1\left(
\frac{
c_1-\Phi^{-1}(0.05)\sqrt{k_1}}{\mu_1}
\right)^2}{
alpha_1>1+beta_1(
(c_1-Phi^{-1}(0.05)sqrt{k_1})/mu_1
\right)^2
}
based on
\deqn{
Pr\left(\frac{\mu_g}{\tau_g^{-1}}\leq c_1 | \tau_g^{-1}\right)<0.05,
}{
Pr(mu_g/tau_g^{-1} <= c_1 | tau_g^{-1})<0.05,
}
where \eqn{c_1=\Phi^{-1}(0.05)}{c_1=Phi^{-1}(0.05)}
and \eqn{\Phi}{Phi} is the cumulative distribution function 
of the standard normal distribution.

For cluster 2, we add one constraint
\deqn{
\alpha_2>1+\beta_2\left(
\frac{
c_2-\Phi^{-1}(0.95)\sqrt{k_2}}{\mu_2}
\right)^2}{
alpha_2>1+beta_2(
(c_2-Phi^{-1}(0.95)sqrt{k_2})/mu_2
\right)^2
}
based on
\deqn{
Pr\left(\frac{\mu_g}{\tau_g^{-1}} \geq c_2 | \tau_g^{-1}\right)<0.05,
}{
Pr(mu_g/tau_g^{-1} >= c_2 | tau_g^{-1})<0.05,
}
where \eqn{c_2=\Phi^{-1}(0.95)}{c_2=Phi^{-1}(0.95)}
and \eqn{\Phi}{Phi} is the cumulative distribution function 
of the standard normal distribution.

To do unconstraint numerical optimization, we do parameter
reparameterization:
\deqn{
\mu_1=\exp(\delta_1),
k_1=\Phi(\xi_1),
\beta_1=\exp(\nu_1),\\
\alpha_1=\exp(\lambda_1)+1+\beta_1\left(
\frac{c_1-\Phi^{-1}(0.05)\sqrt{k_1}}{\mu_1}
\right)^2,\\
\mu_2= -\exp(\delta_2),
k_2=\Phi(\xi_2),
\beta_2=\exp(\nu_2),\\
\alpha_2=\exp(\lambda_2)+1+\beta_2\left(
\frac{c_2-\Phi^{-1}(0.95)\sqrt{k_2}}{\mu_2}
\right)^2,\\
\beta_3=\exp(\nu_3),
\alpha_3=\exp(\lambda_3).
}{
mu_1=exp(delta_1),
k_1=Phi(xi_1),
beta_1=exp(nv_1),\\
alpha_1=exp(lambda_1)+1+beta_1left(
frac{c_1-Phi^{-1}(0.05)sqrt{k_1}}{mu_1}
right)^2,\\
mu_2= -exp(delta_2),
k_2=Phi(xi_2),
beta_2=exp(nv_2),\\
alpha_2=exp(lambda_2)+1+beta_2left(
frac{c_2-Phi^{-1}(0.95)sqrt{k_2}}{mu_2}
right)^2,\\
beta_3=exp(nv_3),
alpha_3=exp(lambda_3).
}

}
\value{

A \eqn{10\times 1}{10x1} vector of reparameterized parameters:
\eqn{\mu_1}{mu_1},
\eqn{k_1}{k_1},
\eqn{\alpha_1}{alpha-1},
\eqn{\beta_1}{beta_1},
\eqn{\alpha_3}{alpha_3},
\eqn{\beta_3}{beta_3},
}
\references{Li Y, Morrow J, Raby B, Tantisira K, Weiss ST, Huang W, Qiu W. (2017), <doi:10.1371/journal.pone.0174602>}
\author{
Yunfeng Li <colinlee1999@gmail.com> and Weiliang Qiu <stwxq@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

\seealso{
 See Also as \code{\link{getRePara}}
}
\examples{
getPara.orig(
  delta1 = -0.690142787, 
  xi1 = -7.212004793, 
  lambda1 = -13.152520780, 
  nu1 = -2.199687707,
  delta2 = -0.168584053, 
  xi2 = 0.008683666, 
  lambda2 = -13.582936416, 
  nu2 = -2.671150369,
  lambda3 = 0.331454152, 
  nu3 = -2.339660241,
  c1 = qnorm(0.95), 
  c2 = qnorm(0.05)
)

}
\keyword{ method }% use one of  RShowDoc("KEYWORDS")

