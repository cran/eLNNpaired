\name{checkPara}
\alias{checkPara}
\title{
Check Parameter Values
}
\description{
Check if parameter values satisfy constraints.
}
\usage{
checkPara(
  mu1, 
  k1, 
  alpha1, 
  beta1, 
  mu2, 
  k2, 
  alpha2, 
  beta2, 
  alpha3, 
  beta3, 
  c1 = qnorm(0.95), 
  c2 = qnorm(0.05))
}
\arguments{
  \item{mu1}{
mean of the mean expression levels for gene probes in cluster 1 (over-expressed probes).
}
  \item{k1}{
scalar in the variance of the mean expression levles for gene probes
in cluster 1 (over-expressed probes).
}
  \item{alpha1}{
shape parameter of the distribution of the variance of gene expression
levels for gene probes in cluster 1 (over-expressed probes).
}
  \item{beta1}{
rate parameter of the distribution of the variance of gene expression
levels for gene probes in cluster 1 (over-expressed probes).
}
  \item{mu2}{
mean of the mean expression levels for gene probes in cluster 2 (under-expressed probes).
}
  \item{k2}{
scalar in the variance of the mean expression levles for gene probes
in cluster 2 (under-expressed probes).
}
  \item{alpha2}{
shape parameter of the distribution of the variance of gene expression
levels for gene probes in cluster 2 (under-expressed probes).
}
  \item{beta2}{
rate parameter of the distribution of the variance of gene expression
levels for gene probes in cluster 2 (under-expressed probes).
}
  \item{alpha3}{
shape parameter of the distribution of the variance of gene expression
levels for gene probes in cluster 3 (non-differentially expressed probes).
}
  \item{beta3}{
rate parameter of the distribution of the variance of gene expression
levels for gene probes in cluster 3 (non-differentially expressed probes).
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
beta_1=exp(nu_1),\\
alpha_1=exp(lambda_1)+1+beta_1left(
frac{c_1-Phi^{-1}(0.05)sqrt{k_1}}{mu_1}
right)^2,\\
mu_2= -exp(delta_2),
k_2=Phi(xi_2),
beta_2=exp(nu_2),\\
alpha_2=exp(lambda_2)+1+beta_2left(
frac{c_2-Phi^{-1}(0.95)sqrt{k_2}}{mu_2}
right)^2,\\
beta_3=exp(nu_3),
alpha_3=exp(lambda_3).
}

}
\value{
Return zero if parameter values satisfy constraints. Otherwise,
program exits abnormally.
}
\references{Li Y, Morrow J, Raby B, Tantisira K, Weiss ST, Huang W, Qiu W. (2017), <doi:10.1371/journal.pone.0174602>}
\author{
Yunfeng Li <colinlee1999@gmail.com> and Weiliang Qiu <stwxq@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{

checkPara(
  mu1 = 5.015045e-01, 
  k1 = 2.756702e-13, 
  alpha1 = 2.192323e+00, 
  beta1 = 1.108378e-01,
  mu2 = -8.448602e-01, 
  k2 = 5.034642e-01, 
  alpha2 = 1.766274e+00, 
  beta2 = 6.917261e-02,
  alpha3 = 1.392992e+00, 
  beta3 = 9.636037e-02,
  c1 = qnorm(0.95),
  c2 = qnorm(0.05)
)


}
\keyword{ method }% use one of  RShowDoc("KEYWORDS")

