# created by Weiliang Qiu on Jan. 3, 2017

# translate re-parameterized parameters to original scale
getPara.orig=function(
  delta1, xi1, lambda1, nu1,
  delta2, xi2, lambda2, nu2,
  lambda3, nu3,
  c1=qnorm(0.95), 
  c2=qnorm(0.05)
)
{
  mu1=exp(delta1)
  k1=pnorm(xi1)
  beta1=exp(nu1)
  numer1=c1-qnorm(0.05)*sqrt(k1)
  alpha1=exp(lambda1)+1+beta1*(numer1/mu1)^2

  mu2= -exp(delta2)
  k2=pnorm(xi2)
  beta2=exp(nu2)
  numer2=c2-qnorm(0.95)*sqrt(k2)
  alpha2=exp(lambda2)+1+beta2*(numer2/mu2)^2

  beta3=exp(nu3)
  alpha3=exp(lambda3)

  res=c(
    mu1, k1, alpha1, beta1,
    mu2, k2, alpha2, beta2,
    alpha3, beta3)
  names(res)=c(
    "mu1", "k1", "alpha1", "beta1",
    "mu2", "k2", "alpha2", "beta2",
    "alpha3", "beta3")

  return(res)
}

# check if parameters satisfy constraints
checkPara=function(
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
  c2 = qnorm(0.05)
)
{
  if(mu1<0)
  {
    stop("mu1 should be positive!\n")
  }
  if(mu2>0)
  {
    stop("mu2 should be negative\n")
  }

  if(k1<0)
  {
    stop("k1 should be positive!\n")
  }
  if(k2<0)
  {
    stop("k2 should be positive!\n")
  }

  if(alpha1<0)
  {
    stop("alpha1 should be positive!\n")
  }
  if(beta1<0)
  {
    stop("beta1 should be positive!\n")
  }

  if(alpha2<0)
  {
    stop("alpha2 should be positive!\n")
  }
  if(beta2<0)
  {
    stop("beta2 should be positive!\n")
  }

  if(alpha3<0)
  {
    stop("alpha3 should be positive!\n")
  }
  if(beta3<0)
  {
    stop("beta3 should be positive!\n")
  }

  numer1=c1 - qnorm(0.05)*sqrt(k1)
  cut1=1+beta1*(numer1/mu1)^2

  if(alpha1 <= cut1)
  {
    stop("alpha1 must be >", cut1, "\n")
  }

  numer2=c2 - qnorm(0.95)*sqrt(k2)
  cut2=1+beta2*(numer2/mu2)^2

  if(alpha2 <= cut2)
  {
    stop("alpha2 must be >", cut2, "\n")
  }

  return(0)
}


# re-parameterized parameters
getRePara=function(
  mu1, k1, alpha1, beta1,
  mu2, k2, alpha2, beta2,
  alpha3, beta3,
  c1=qnorm(0.95), 
  c2=qnorm(0.05)
)
{
  delta1=log(mu1)
  xi1=qnorm(k1)
  nu1=log(beta1)
  numer1=c1-qnorm(0.05)*sqrt(k1)
  lambda1=log(alpha1-1-beta1*(numer1/mu1)^2)

  delta2=log(-mu2)
  xi2=qnorm(k2)
  nu2=log(beta2)
  numer2=c2-qnorm(0.95)*sqrt(k2)
  lambda2=log(alpha2-1-beta2*(numer2/mu2)^2)

  nu3=log(beta3)
  lambda3=log(alpha3)

  res=c(
    delta1, xi1, lambda1, nu1,
    delta2, xi2, lambda2, nu2,
    lambda3, nu3)
  names(res)=c(
    "delta1", "xi1", "lambda1", "nu1",
    "delta2", "xi2", "lambda2", "nu2",
    "lambda3", "nu3")

  return(res)
}

