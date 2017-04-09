#
# 1
# #' Calculate log f1, log f2 and log f3.
# #'
# #' @param psi It contains the paramters for each likelihood function.
# #' @param sum_dgl_by_l An G * 1 matrix, the summation result of every row of 'exprs(E_Set)'.
# #' @param sum_dgl_square_by_l An G * 1 matrix, the summation of every squared elements of 'exprs(E_Set)' by row.
# #' @param n The number of samples.
# #' @param c1 Parameter used in constraints.
# #' @param c2 Parameter used in constraints.
# #' @return A G by 3 matrix, first column for log f1, second column for log f2 and third column for log f3.
lf123 <- function(
  psi,
  sum_dgl_by_l,
  sum_dgl_square_by_l,
  n,
  c1,
  c2)
{
  #result = matrix(rep(0,length(sum_dgl_by_l)*3),length(sum_dgl_by_l),3)
  result = matrix(0, nrow=length(sum_dgl_by_l), ncol=3)
  colnames(result) = c("logf1", "logf2", "logf3")

#  para.orig=getPara.orig(
#    delta1=psi[1], xi1=psi[2], lambda1=psi[3], nu1=psi[4],
#    delta2=psi[5], xi2=psi[6], lambda2=psi[7], nu2=psi[8],
#    lambda3=psi[9], nu3=psi[10],
#    c1=c1, 
#    c2=c2
#  )
#

  ####
  # cluster 1
  ####
  mu_0 = exp(psi[1])
  k = pnorm(psi[2])
  beta = exp(psi[4])
  alpha = exp(psi[3]) + 1 + beta * ((c1 - qnorm(0.05) * sqrt(k))/mu_0)^2

#  mu_0 = para.orig[1]
#  k = para.orig[2]
#  alpha = para.orig[3]
#  beta = para.orig[4]
#
  A = n/(2*(n*k+1))
  B_g = sum_dgl_by_l/n
  C_g = beta + sum_dgl_square_by_l/2 - (sum_dgl_by_l)^2/(2*n)
  D = alpha * log(beta) + lgamma(n/2+alpha) - lgamma(alpha) - n/2 * log(2*pi) - log(n*k+1)/2
  result[,1] = D - (n/2+alpha) * log (C_g + A * (mu_0 - B_g)^2)

  ####
  # cluster 2
  ####
  mu_0 = -exp(psi[5])
  k = pnorm(psi[6])
  beta = exp(psi[8])
  alpha = exp(psi[7]) + 1 + beta * ((c2 - qnorm(0.95) * sqrt(k))/mu_0)^2

#  mu_0 = para.orig[5]
#  k = para.orig[6]
#  alpha = para.orig[7]
#  beta = para.orig[8]
#

  A = n/(2*(n*k+1))
  B_g = sum_dgl_by_l/n
  C_g = beta + sum_dgl_square_by_l/2 - (sum_dgl_by_l)^2/(2*n)
  D = alpha * log(beta) + lgamma(n/2+alpha) - lgamma(alpha) - n/2 * log(2*pi) - log(n*k+1)/2
  result[,2] = D - (n/2+alpha) * log (C_g + A * (mu_0 - B_g)^2)

  ####
  # cluster 3
  ####
  beta = exp(psi[10])
  alpha = exp(psi[9])

  #alpha = para.orig[9]
  #beta = para.orig[10]

  D = alpha * log(beta) + lgamma(n/2+alpha) - lgamma(alpha) - n/2 * log(2*pi)
  result[,3] = D - (n/2+alpha) * log (beta + sum_dgl_square_by_l/2)

  return (result)
}

# 2
# #' Calculate tilde z.
# #'
# #' @param psi It contains the paramters for each likelihood function.
# #' @param t_pi, It is 'pi' as in paper, we do not use 'pi' since it is reserved as a mathematical constatnt in R.
# #' @param sum_dgl_by_l An G * 1 matrix, the summation result of every row of 'exprs(E_Set)'.
# #' @param sum_dgl_square_by_l An G * 1 matrix, the summation of every squared elements of 'exprs(E_Set)' by row.
# #' @param n The number of samples.
# #' @param c1 Parameter used in constraints.
# #' @param c2 Parameter used in constraints.
# #' @return A G by 3 matrix.
get_tilde_z <- function(
  psi,
  t_pi,
  sum_dgl_by_l,
  sum_dgl_square_by_l,
  n,
  c1,
  c2)
{
  logf = lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, n, c1, c2)
  max_logf = apply(logf, 1, max, na.rm = TRUE)

  pi1=t_pi[1]
  pi2=t_pi[2]
  pi3=1-pi1-pi2

  t1 = pi1 * exp(logf[,1] - max_logf)
  t2 = pi2 * exp(logf[,2] - max_logf)
  t3 = pi3 * exp(logf[,3] - max_logf)
  total = t1 + t2 + t3

  result = cbind(t1, t2, t3)/total

  return(result)
}

# 3
# #' Calculate negative value of the expected complete data log likelihood.
# #'
# #' @param psi It contains the paramters for each likelihood function.
# #' @param t_pi, It is 'pi' as in paper, we do not use 'pi' since it is reserved as a mathematical constatnt in R.
# #' @param sum_dgl_by_l An G * 1 matrix, the summation result of every row of 'exprs(E_Set)'.
# #' @param sum_dgl_square_by_l An G * 1 matrix, the summation of every squared elements of 'exprs(E_Set)' by row.
# #' @param n The number of samples.
# #' @param tilde_z Tilde z as in paper.
# #' @param b Parameter for Dirichlet distribution.
# #' @param c1 Parameter used in constraints.
# #' @param c2 Parameter used in constraints.
# #' @return A single variable indicating the value of complete data log likelihood.
#l_c <- function(
negative_l_c <- function(
  psi,
  t_pi,
  sum_dgl_by_l,
  sum_dgl_square_by_l,
  n,
  tilde_z,
  b,
  c1,
  c2)
{
  logf = lf123(psi, sum_dgl_by_l, sum_dgl_square_by_l, n, c1, c2)

  pi1=t_pi[1]
  pi2=t_pi[2]
  pi3=1-pi1-pi2

  result = 0
  result = result + sum(tilde_z[,1] * logf[,1], na.rm = TRUE)
  result = result + sum(tilde_z[,2] * logf[,2], na.rm = TRUE)
  result = result + sum(tilde_z[,3] * logf[,3], na.rm = TRUE)
  result = result + sum(tilde_z[,1] * log(pi1), na.rm = TRUE)
  result = result + sum(tilde_z[,2] * log(pi2), na.rm = TRUE)
  result = result + sum(tilde_z[,3] * log(pi3), na.rm = TRUE)

  result = result + lgamma(b[1]+b[2]+b[3]) - lgamma(b[1]) - lgamma(b[2]) - lgamma(b[3])
  t_pi2=c(pi1, pi2, pi3)
  result = result + sum((b-1) * log(t_pi2), na.rm = TRUE)

  return ( -result)
}

## 4
## #' Calculate negative expected complete data log likelihood.
## #'
## #' @param psi It contains the paramters for each likelihood function.
## #' @param t_pi, It is 'pi' as in paper, we do not use 'pi' since it is reserved as a mathematical constatnt in R.
## #' @param sum_dgl_by_l An G * 1 matrix, the summation result of every row of 'exprs(E_Set)'.
## #' @param sum_dgl_square_by_l An G * 1 matrix, the summation of every squared elements of 'exprs(E_Set)' by row.
## #' @param n The number of samples.
## #' @param tilde_z Tilde z as in paper.
## #' @param b Parameter for Dirichlet distribution.
## #' @param c1 Parameter used in constraints.
## #' @param c2 Parameter used in constraints.
## #' @return A single variable indicating the negative value of complete data log likelihood.
#negative_l_c <- function(
#  psi,
#  t_pi,
#  sum_dgl_by_l,
#  sum_dgl_square_by_l,
#  n,
#  tilde_z,
#  b,
#  c1,
#  c2)
#{
#  return (-l_c(
#    psi,
#    t_pi,
#    sum_dgl_by_l,
#    sum_dgl_square_by_l,
#    n,
#    tilde_z,
#    b,
#    c1,
#    c2))
#}

# 5
# #' Calculate gradient of l_c w.r.t. psi.
# #'
# #' @param psi It contains the paramters for each likelihood function.
# #' @param t_pi, It is 'pi' as in paper, we do not use 'pi' since it is reserved as a mathematical constatnt in R.
# #' @param sum_dgl_by_l An G * 1 matrix, the summation result of every row of 'exprs(E_Set)'.
# #' @param sum_dgl_square_by_l An G * 1 matrix, the summation of every squared elements of 'exprs(E_Set)' by row.
# #' @param n The number of samples.
# #' @param tilde_z Tilde z as in paper.
# #' @param b Parameter for Dirichlet distribution.
# #' @param c1 Parameter used in constraints.
# #' @param c2 Parameter used in constraints.
# #' @return A vector that is the same length as psi, containing the gradient.
#gradient_l_c <- function(
gradient_negative_l_c <- function(
  psi,
  t_pi,
  sum_dgl_by_l,
  sum_dgl_square_by_l,
  n,
  tilde_z,
  b,
  c1,
  c2)
{
  G = length(sum_dgl_by_l)
  #gradient_temp = matrix(rep(0,12),4,3)

#  para.orig=getPara.orig(
#    delta1=psi[1], xi1=psi[2], lambda1=psi[3], nu1=psi[4],
#    delta2=psi[5], xi2=psi[6], lambda2=psi[7], nu2=psi[8],
#    lambda3=psi[9], nu3=psi[10],
#    c1=c1, 
#    c2=c2
#  )
#

  ####
  # cluster 1
  ####
  cluster = 1

  mu_0 = exp(psi[1])
  k = pnorm(psi[2])
  beta = exp(psi[4])
  lambda = psi[3]
  alpha = exp(lambda) + 1 + beta * ((c1 - qnorm(0.05) * sqrt(k))/mu_0)^2

#  mu_0 = para.orig[1]
#  k = para.orig[2]
#  alpha = para.orig[3]
#  beta = para.orig[4]
#  lambda = psi[3]
#
  A = n/(2*(n*k+1))
  B_g = sum_dgl_by_l/n
  C_g = beta + sum_dgl_square_by_l/2 - (sum_dgl_by_l)^2/(2*n)
  D = alpha * log(beta) + lgamma(n/2+alpha) - lgamma(alpha) -n/2 * log(2*pi) - log(n*k+1)/2

  d_mu_0 = - sum(tilde_z[,cluster] * (2*A*(mu_0 - B_g))/(A*(mu_0 - B_g)^2 + C_g) * (n/2 + alpha))
  d_k = - (n * sum(tilde_z[,cluster])) / (2*(n*k+1)) + (n^2)/(2*(n*k+1)^2) * (n/2+alpha) *
          sum(tilde_z[,cluster]*(mu_0-B_g)^2/(A*(mu_0-B_g)^2+C_g))
  d_beta = sum(tilde_z[,cluster]) * alpha / beta -  (n/2 + alpha) *
          (sum(tilde_z[,cluster] / (A * (mu_0 - B_g)^2 + C_g)))
  d_alpha = sum(tilde_z[,cluster] * (log(beta) + digamma(n/2+alpha) - digamma(alpha))) -
          (sum(tilde_z[,cluster] * log(A * (mu_0 - B_g)^2 + C_g)))

  d_delta = mu_0 * (d_mu_0 - 2 * d_alpha * (beta*c1^2*(1+sqrt(k))^2)/mu_0^3)
  d_nu = beta * (d_beta + d_alpha * c1^2 * (1+sqrt(k))^2/mu_0^2)
  d_xi = dnorm(qnorm(k)) * (d_k + d_alpha * beta * c1^2 * (1 + 1/sqrt(k))/mu_0^2)
  d_lambda = exp(lambda) * d_alpha

  gradient_temp = c(d_delta, d_xi, d_lambda, d_nu)

  ####
  # cluster 2
  ####
  cluster = 2

  mu_0 = -exp(psi[5])
  k = pnorm(psi[6])
  beta = exp(psi[8])
  lambda = psi[7]
  alpha = exp(lambda) + 1 + beta * ((c2 - qnorm(0.95) * sqrt(k))/mu_0)^2

#  mu_0 = para.orig[5]
#  k = para.orig[6]
#  alpha = para.orig[7]
#  beta = para.orig[8]
#  lambda = psi[7]

  A = n/(2*(n*k+1))
  B_g = sum_dgl_by_l/n
  C_g = beta + sum_dgl_square_by_l/2 - (sum_dgl_by_l)^2/(2*n)
  D = alpha * log(beta) + lgamma(n/2+alpha) - lgamma(alpha) -n/2 * log(2*pi) - log(n*k+1)/2

  d_mu_0 = - sum(tilde_z[,cluster] * (2*A*(mu_0 - B_g))/(A*(mu_0 - B_g)^2 + C_g) * (n/2 + alpha))
  d_k = - (n * sum(tilde_z[,cluster])) / (2*(n*k+1)) + (n^2)/(2*(n*k+1)^2) * (n/2+alpha) *
          sum(tilde_z[,cluster]*(mu_0-B_g)^2/(A*(mu_0-B_g)^2+C_g))
  d_beta = sum(tilde_z[,cluster]) * alpha / beta -  (n/2 + alpha) *
          (sum(tilde_z[,cluster] / (A * (mu_0 - B_g)^2 + C_g)))
  d_alpha = sum(tilde_z[,cluster] * (log(beta) + digamma(n/2+alpha) - digamma(alpha))) -
          (sum(tilde_z[,cluster] * log(A * (mu_0 - B_g)^2 + C_g)))

  d_delta = mu_0 * (d_mu_0 - 2 * d_alpha * (beta*c2^2*(1+sqrt(k))^2)/mu_0^3)
  d_nu = beta * (d_beta + d_alpha * c2^2 * (1+sqrt(k))^2/mu_0^2)
  d_xi = dnorm(qnorm(k)) * (d_k + d_alpha * beta * c2^2 * (1 + 1/sqrt(k))/mu_0^2)
  d_lambda = exp(lambda) * d_alpha

  gradient_temp = c(gradient_temp, d_delta, d_xi, d_lambda, d_nu)

  ####
  # cluster 3
  ####
  cluster = 3

  mu_0 = 0
  beta = exp(psi[10])
  lambda = psi[9]
  alpha = exp(lambda)

#  mu_0 = 0
#  alpha = para.orig[9]
#  beta = para.orig[10]
#  lambda = psi[9]
#

  d_beta = sum(tilde_z[,cluster]) * alpha / beta -  (n/2 + alpha) *
          (sum(tilde_z[,cluster] / (beta + sum_dgl_square_by_l/2)))
  d_alpha = sum(tilde_z[,cluster] * (log(beta) + digamma(n/2+alpha) - digamma(alpha))) -
          (sum(tilde_z[,cluster] * log(beta + sum_dgl_square_by_l/2)))

  #d_delta = 0
  #d_xi = 0
  d_nu = beta * d_beta
  d_lambda = exp(lambda) * d_alpha

  #gradient_temp[,3] = c(d_delta, d_xi, d_lambda, d_nu)
  gradient_temp = c(gradient_temp, d_lambda, d_nu)

  result = gradient_temp

  return ( -result)
}

## 6
## #' Calculate negative gradient of l_c w.r.t. psi.
## #'
## #' @param psi It contains the paramters for each likelihood function.
## #' @param t_pi, It is 'pi' as in paper, we do not use 'pi' since it is reserved as a mathematical constatnt in R.
## #' @param sum_dgl_by_l An G * 1 matrix, the summation result of every row of 'exprs(E_Set)'.
## #' @param sum_dgl_square_by_l An G * 1 matrix, the summation of every squared elements of 'exprs(E_Set)' by row.
## #' @param n The number of samples.
## #' @param tilde_z Tilde z as in paper.
## #' @param b Parameter for Dirichlet distribution.
## #' @param c1 Parameter used in constraints.
## #' @param c2 Parameter used in constraints.
## #' @return A vector that is the same length as psi, containing the negative gradient.
#gradient_negative_l_c <- function(
#  psi,
#  t_pi,
#  sum_dgl_by_l,
#  sum_dgl_square_by_l,
#  n,
#  tilde_z,
#  b,
#  c1,
#  c2)
#{
#  return (-gradient_l_c(
#    psi,
#    t_pi,
#    sum_dgl_by_l,
#    sum_dgl_square_by_l,
#    n,
#    tilde_z,
#    b,
#    c1,
#    c2))
#}
#
# 7
# #' Calculate pi as in paper.
# #'
# #' @param psi It contains the paramters for each likelihood function.
# #' @param t_pi, It is 'pi' as in paper, we do not use 'pi' since it is reserved as a mathematical constatnt in R.
# #' @param sum_dgl_by_l An G * 1 matrix, the summation result of every row of 'exprs(E_Set)'.
# #' @param sum_dgl_square_by_l An G * 1 matrix, the summation of every squared elements of 'exprs(E_Set)' by row.
# #' @param n The number of samples.
# #' @param tilde_z Tilde z as in paper.
# #' @param b Parameter for Dirichlet distribution.
# #' @return A G by 3 matrix.
get_t_pi <- function(
  psi,
  t_pi,
  sum_dgl_by_l,
  sum_dgl_square_by_l,
  n,
  tilde_z,
  b)
{
  denominator = length(sum_dgl_by_l) + sum(b, na.rm=TRUE) - 3

  t1 = (sum(tilde_z[,1], na.rm=TRUE) + b[1] - 1) / denominator
  t2 = (sum(tilde_z[,2], na.rm=TRUE) + b[2] - 1) / denominator
  t3 = (sum(tilde_z[,3], na.rm=TRUE) + b[3] - 1) / denominator

  return (c(t1,t2,t3))
}

# 8
#' Implementation of eLNNpaired model.
#'
#' @return A list of 5 elementes: E_Set, an ExpressionSet, added clustering results to fData(E_Set)$memGenes in which E_Set is the original input; mleinfo, the return value from built-in optim function; t_pi, the estimation of 'pi' as in paper; repeated_times, the number of iterations that EM algorithm takes; tilde_z, tilde z as discribed in paper, the probability that each gene will fall into each cluster.
#' @param E_Set An ExpressionSet on which model-based clustering will be performed.
#' @param b A vector of concentration parameters used in Dirichlet distribution. Default value is b = c(2,2,2).
#' @param verbose An indicator variable telling if print out intermediate results: zero value for not printing out, non-zero for printing out. Default value is verbose = 0.
#' @param converge_threshold One of the two termination criteria of iteration. The smaller this value is set, the harder the optimization procedure in eLNNpaired will be considered to be converged. Default value is converge_threshold = 1e-6.
#' @param using_limit An indicator variable telling if constrained optimization method will be adopted: zero value for unconstrained while non-zero value for constrained optimization. Default value is using_limit = 0. Note that even if using_limit is set to be 0, the following two parameters param_limit_min and param_limit_max may still be used when initial estimation of one or more of the parameters are not feasible.
#' @param param_limit_min An vector of lower bounds of parameters. Default value is param_limit_min = c(-6,qnorm(0.01),-6,-6,-6,qnorm(0.01),-6,-6,-6,-6).
#' @param param_limit_max An vector of upper bounds of parameters. Default value is param_limit_max = c(6,qnorm(0.99),6,6,6,qnorm(0.99),6,6,6,6) = -param_limit_min.
#' @param maxIT An integer, the max allowed number of iterations in R built-in function optim. Default value is maxIT = 100.
#' @param maxRT An integer, the max allowed number of iterations for EM algorithm. Default value is maxRT = 500.
#' @param c1 A parameter in constraints. It should be in the form of c1 = qnorm(X), where X is a decimal smaller than 1 but close to 1. Larger X gives more stringent constraint. Default value is c1 = qnorm(0.95).
#' @param c2 A parameter in constraints. It should be in the form of c2 = qnorm(Y), where Y is a decimal larger than 0 but close to 0. Smaller Y gives more stringent constraint. Default value is c2 = qnorm(0.05).
#' @export
#' @example example/example.r
eLNNpaired <- function(
  E_Set,
  b = c(2,2,2),
  verbose = 0,
  converge_threshold = 1e-6,
  using_limit = 0,
  param_limit_min = c(-6,qnorm(0.01),-6,-6,-6,qnorm(0.01),-6,-6,-6,-6),
  param_limit_max = c(6,qnorm(0.99),6,6,6,qnorm(0.99),6,6,6,6),
  maxIT = 100,
  maxRT = 100,
  c1 = qnorm(0.95),
  c2 = qnorm(0.05)
  )
{
  # 'G' is the number of genes
  # 'n' is the number of samples
  G = nrow(E_Set)
  n = ncol(E_Set)

  data_matrix_of_E_Set = exprs(E_Set)

  # 'sum_dgl_by_l' is an G * 1 matrix, the summation result of every row of 'data_matrix_of_E_Set'
  sum_dgl_by_l = apply(data_matrix_of_E_Set,1,sum, na.rm=TRUE)

  # 'sum_dgl_square_by_l' is an G * 1 matrix, the summation of every squared elements of 'data_matrix_of_E_Set' by row
  sum_dgl_square_by_l = apply(data_matrix_of_E_Set^2,1,sum, na.rm = TRUE)

  #column_names = colnames(fData(E_Set))
  nprobes=nrow(E_Set)
  probeid=1:nprobes
  genes=1:nprobes
  chr=rep(1, nprobes)
  fData(E_Set)$myprobeid=probeid
  fData(E_Set)$mygenes=genes
  fData(E_Set)$mychr=chr

  result_limma = lmFitPaired(
    E_Set,
    probeID.var = "myprobeid",
    gene.var = "mygenes",
    chr.var = "mychr",
    verbose = FALSE)

  frame_unsorted = result_limma$frame.unsorted
  # over_expressed_sub_script = frame_unsorted$pos[which(frame_unsorted$stats > 0 & frame_unsorted$p.adj < 0.05)]
  # under_expressed_sub_script = frame_unsorted$pos[which(frame_unsorted$stats < 0 & frame_unsorted$p.adj < 0.05)]
  # non_diff_sub_script = frame_unsorted$pos[which(frame_unsorted$p.adj >= 0.05)]

  over_expressed_sub_script = frame_unsorted$pos[which(frame_unsorted$stats > 0 & frame_unsorted$pval < 0.05)]
  under_expressed_sub_script = frame_unsorted$pos[which(frame_unsorted$stats < 0 & frame_unsorted$pval < 0.05)]
  non_diff_sub_script = frame_unsorted$pos[which(frame_unsorted$pval >= 0.05)]

  pi_prior = c(
    length(under_expressed_sub_script)/G,
    length(over_expressed_sub_script)/G,
    0)
  pi_prior[3] = 1 - pi_prior[1] - pi_prior[2]

  over_expressed_E_Set = E_Set[over_expressed_sub_script,]
  under_expressed_E_Set = E_Set[under_expressed_sub_script,]
  non_diff_E_Set = E_Set[non_diff_sub_script,]

  #########################################
  # cluster 1
  data_matrix_of_E_Set = exprs(over_expressed_E_Set)

  median_dgl_by_l = apply(data_matrix_of_E_Set, 1, median, na.rm=TRUE)

  sorted_median_dgl_by_l = sort(median_dgl_by_l)

  temp = median(sorted_median_dgl_by_l, na.rm=TRUE)
  if (temp>0) delta_1 = log(temp)
  else delta_1 = param_limit_min[1]

  temp_tau = 1 / (apply(data_matrix_of_E_Set, 1, mad, na.rm=TRUE)^2)
  omega = mad(median_dgl_by_l)^2
  k_prior = omega * median(temp_tau, na.rm=TRUE)

  if (k_prior>=1) xi_1 = param_limit_max[2]
  else xi_1 = qnorm(k_prior)

  alpha = (median(temp_tau, na.rm=TRUE)^2)/(mad(temp_tau, na.rm=TRUE)^2)
  beta = median(temp_tau, na.rm=TRUE)/(mad(temp_tau, na.rm=TRUE)^2)
  nu_1 = log(beta)
  exp_lambda_1 = alpha - 1 - beta*c1^2*(1+sqrt(k_prior))^2/exp(delta_1)^2
  if (exp_lambda_1<0) lambda_1 = param_limit_min[3]
  else lambda_1 = log(exp_lambda_1)

  # end of cluster 1
  #########################################

  #########################################
  # cluster 2
  data_matrix_of_E_Set = exprs(under_expressed_E_Set)

  median_dgl_by_l = apply(data_matrix_of_E_Set, 1, median, na.rm=TRUE)

  sorted_median_dgl_by_l = sort(median_dgl_by_l)

  temp = median(sorted_median_dgl_by_l, na.rm=TRUE)
  temp = - temp
  if (temp>0) delta_2 = log(temp)
  else delta_2 = param_limit_min[5]

  temp_tau = 1 / (apply(data_matrix_of_E_Set, 1, mad, na.rm=TRUE)^2)
  omega = mad(median_dgl_by_l)^2
  k_prior = omega * median(temp_tau, na.rm=TRUE)

  if (k_prior>=1) xi_2 = param_limit_max[6]
  else xi_2 = qnorm(k_prior)

  alpha = (median(temp_tau, na.rm=TRUE)^2)/(mad(temp_tau, na.rm=TRUE)^2)
  beta = median(temp_tau, na.rm=TRUE)/(mad(temp_tau, na.rm=TRUE)^2)
  nu_2 = log(beta)
  exp_lambda_2 = alpha - 1 - beta*c2^2*(1+sqrt(k_prior))^2/exp(delta_2)^2
  if (exp_lambda_2<0) lambda_2 = param_limit_min[7]
  else lambda_2 = log(exp_lambda_2)
  # end of cluster 2
  #########################################

  #########################################
  # cluster 3
  data_matrix_of_E_Set = exprs(non_diff_E_Set)

  median_dgl_by_l = apply(data_matrix_of_E_Set, 1, median, na.rm=TRUE)

  sorted_median_dgl_by_l = sort(median_dgl_by_l)

  temp_tau = 1 / (apply(data_matrix_of_E_Set, 1, mad, na.rm=TRUE)^2)

  lambda_3 = log((median(temp_tau, na.rm=TRUE)^2)/(mad(temp_tau, na.rm=TRUE)^2))
  nu_3 = log(median(temp_tau, na.rm=TRUE)/(mad(temp_tau, na.rm=TRUE)^2))
  # end of cluster 3
  #########################################

  data_matrix_of_E_Set = exprs(E_Set)

  psi = c(delta_1, xi_1, lambda_1, nu_1,
          delta_2, xi_2, lambda_2, nu_2,
          lambda_3, nu_3)

  if (verbose)
  {
    print(psi)
    print(pi_prior)
  }

  #t_pi = pi_prior
  t_pi = pi_prior[1:2]

  tilde_z = get_tilde_z(
    psi,
    t_pi,
    sum_dgl_by_l,
    sum_dgl_square_by_l,
    n,
    c1,
    c2)

  if (using_limit)
  {
    mleinfo = optim(par = psi, fn = negative_l_c, gr = gradient_negative_l_c,
      t_pi = t_pi, sum_dgl_by_l = sum_dgl_by_l, 
      sum_dgl_square_by_l = sum_dgl_square_by_l,
      n = n, tilde_z = tilde_z, method = 'L-BFGS-B',
      b = b, c1 = c1, c2 = c2,
      lower=param_limit_min,
      upper=param_limit_max,
      control = list(maxit = maxIT))
  }
  else
  {
    mleinfo = optim(par = psi, fn = negative_l_c, gr = gradient_negative_l_c,
      t_pi = t_pi, sum_dgl_by_l = sum_dgl_by_l, 
      sum_dgl_square_by_l = sum_dgl_square_by_l,
      n = n, tilde_z = tilde_z, method = 'L-BFGS-B',
      b = b, c1 = c1, c2 = c2,
      control = list(maxit = maxIT))
  }

  t_pi = get_t_pi(
    psi,
    t_pi,
    sum_dgl_by_l,
    sum_dgl_square_by_l,
    n,
    tilde_z,
    b)
  t_pi=t_pi[1:2]

  repeated_times = 0

  if(using_limit)
  {
    while (repeated_times<maxRT)
    {
      repeated_times = repeated_times + 1
      if (verbose)
      {
  
      }
  
      psi = mleinfo$par
  
      tilde_z = get_tilde_z(
        psi,
        t_pi,
        sum_dgl_by_l,
        sum_dgl_square_by_l,
        n,
        c1,
        c2)
  
      if (verbose)
      {
        print(c("repeated times:", repeated_times))
        cat('psi >>\n',psi,'\n')
        cat('pi >>\n',t_pi,'\n')
      }
  
      last_mleinfo = mleinfo
  
      mleinfo = optim(par = psi, fn = negative_l_c, gr = gradient_negative_l_c,
        t_pi = t_pi, sum_dgl_by_l = sum_dgl_by_l, 
        sum_dgl_square_by_l = sum_dgl_square_by_l,
        n = n, tilde_z = tilde_z, method = 'L-BFGS-B',
        b = b, c1 = c1, c2 = c2,
        lower=param_limit_min,
        upper=param_limit_max,
        control = list(maxit = maxIT))
     
      t_pi = get_t_pi(
        psi,
        t_pi,
        sum_dgl_by_l,
        sum_dgl_square_by_l,
        n,
        tilde_z,
        b)
      t_pi = t_pi[1:2]
  
      if (sum(abs(last_mleinfo$par - mleinfo$par))<converge_threshold) break
    }
  
  } 
  else {  
    while (repeated_times<maxRT)
    {
      repeated_times = repeated_times + 1
      if (verbose)
      {
  
      }
  
      psi = mleinfo$par
  
      tilde_z = get_tilde_z(
        psi,
        t_pi,
        sum_dgl_by_l,
        sum_dgl_square_by_l,
        n,
        c1,
        c2)
  
      if (verbose)
      {
        print(c("repeated times:", repeated_times))
        cat('psi >>\n',psi,'\n')
        cat('pi >>\n',t_pi,'\n')
      }
  
      last_mleinfo = mleinfo
  
      mleinfo = optim(par = psi, fn = negative_l_c, gr = gradient_negative_l_c,
        t_pi = t_pi, sum_dgl_by_l = sum_dgl_by_l, 
        sum_dgl_square_by_l = sum_dgl_square_by_l,
        n = n, tilde_z = tilde_z, method = 'L-BFGS-B',
        b = b, c1 = c1, c2 = c2,
        control = list(maxit = maxIT))
  
      t_pi = get_t_pi(
        psi,
        t_pi,
        sum_dgl_by_l,
        sum_dgl_square_by_l,
        n,
        tilde_z,
        b)
      t_pi = t_pi[1:2]
  
      if (sum(abs(last_mleinfo$par - mleinfo$par))<converge_threshold) break
    }
  }

  tilde_z = get_tilde_z(
    psi,
    t_pi,
    sum_dgl_by_l,
    sum_dgl_square_by_l,
    n,
    c1,
    c2)

  memGenes = apply(tilde_z,1,which.max)

  psi= mleinfo$par
  names(psi)=c(
    "delta1", "xi1", "lambda1", "nu1",
    "delta2", "xi2", "lambda2", "nu2",
    "lambda3", "nu3")

  para.orig = getPara.orig(
    delta1=psi[1], 
    xi1=psi[2], 
    lambda1=psi[3], 
    nu1=psi[4],
    delta2=psi[5], 
    xi2=psi[6], 
    lambda2=psi[7], 
    nu2=psi[8],
    lambda3=psi[9], 
    nu3=psi[10],
    c1=c1, 
    c2=c2
  )
  names(para.orig)=c(
    "mu1", "k1", "alpha1", "beta1",
    "mu2", "k2", "alpha2", "beta2",
    "alpha3", "beta3")

  memGenes2=rep(1, length(memGenes))
  pos=which(memGenes==3)
  if(length(pos))
  {
    memGenes2[pos]=0
  }

  result = list(
    mleinfo = mleinfo,
    psi = psi,
    para.orig= para.orig,
    memGenes=memGenes,
    memGenes2=memGenes2,
    t_pi = t_pi,
    repeated_times = repeated_times,
    tilde_z = tilde_z)

  invisible(result)
}

# 9
#' Draw probability density of three clusters.
#'
#' @return No return value.
#' @param psi A vector of length 10. It contains the parameters after reparameterization as illustrated in paper. It has no default value and has to be provided by user.
#' @param t_pi A vector of length 3 and summation equal to 1. It is 'pi' as in paper, we do not use 'pi' since it is reserved as a mathematical constatnt in R. It should be provided when weighted = 1. Default value is t_pi = c(0.3333, 0.3333, 0.3334).
#' @param x A range on which the graph will be drawn. Default value is x = seq(from = -6, to = 6, by = 0.1).
#' @param c1 A parameter in constraints. It should be in the form of c1 = qnorm(X), where X is a decimal smaller than 1 but close to 1. Larger X gives more stringent constraint. Default value is c1 = qnorm(0.95).
#' @param c2 A parameter in constraints. It should be in the form of c2 = qnorm(Y), where Y is a decimal larger than 0 but close to 0. Smaller Y gives more stringent constraint. Default value is c2 = qnorm(0.05).
#' @param weighted An indicator variable telling if draw weighted densities: zero value for drawing unweighted densities, non-zero value for drawing weighted densities. The weights are provided by t_pi. Default value is weighted = 1.
#' @param draw_hist An indicator variable telling if draw histogram from E_Set: zero value for not drawing histogram, non-zero value for drawing histogram. Default value is draw_hist = 0.
#' @param E_Set An ExpressionSet based on which a histgram will be drawn. It will only be used when draw_hist is non-zero. Default value is E_Set = NULL.
#' @export
#' @example example/example.r
draw_density <- function(
  psi,
  t_pi = c(0.3333,0.3333),
  x = seq(from = -6, to = 6, by = 0.1),
  c1 = qnorm(0.95),
  c2 = qnorm(0.05),
  weighted = 1,
  draw_hist = 0,
  E_Set = NULL)
{
  logf = lf123(
          psi,
          x,
          x^2,
          1,
          c1,
          c2)

  pi1=t_pi[1]
  pi2=t_pi[2]
  pi3=1-pi1-pi2

  line_width = 5

  if (draw_hist)
  {
    if (weighted)
    {
      y1 = exp(logf[,1]) * (pi1)
      y2 = exp(logf[,2]) * (pi2)
      y3 = exp(logf[,3]) * (pi3)
      y4 = hist(c(exprs(E_Set)),
        plot = FALSE)
      range1 = range(y1)
      range2 = range(y2)
      range3 = range(y3)
      range4 = range(y4$density)
      plot(y4,
        ylim = c(0, max(range1+range2+range3,range4)),
        main = "weighted probability density of three clusters",
        xlab = "x",
        ylab = "weighted density",
        freq = FALSE)
    }
    else
    {
      y1 = exp(logf[,1])
      y2 = exp(logf[,2])
      y3 = exp(logf[,3])
      y4 = hist(c(exprs(E_Set)),
         plot = FALSE)
      range1 = range(y1)
      range2 = range(y2)
      range3 = range(y3)
      range4 = range(y4$density)
      plot(y4,
        ylim = c(0, max(range1+range2+range3,range4)),
        main = "probability density of three clusters",
        xlab = "x",
        ylab = "density",
        freq = FALSE)
    }

    lines(x = x, y = y1, col = 'red', lwd = line_width)
    lines(x = x, y = y2, col = 'blue', lwd = line_width)
    lines(x = x, y = y3, col = 'black', lwd = line_width)
    lines(x = x, y = y1+y2+y3, col = 'brown', lwd = line_width)
  }
  else
  {
    if (weighted)
    {
      y1 = exp(logf[,1]) * (pi1)
      y2 = exp(logf[,2]) * (pi2)
      y3 = exp(logf[,3]) * (pi3)
      range1 = range(y1)
      range2 = range(y2)
      range3 = range(y3)
      plot(x = x, y = y1,
        ylim = c(0,max(range1,range2,range3)),
        col = 'red',
        lwd = line_width,
        type = 'l',
        main = "weighted probability density of three clusters",
        xlab = "x",
        ylab = "weighted density")
    }
    else
    {
      y1 = exp(logf[,1])
      y2 = exp(logf[,2])
      y3 = exp(logf[,3])
      range1 = range(y1)
      range2 = range(y2)
      range3 = range(y3)
      plot(x = x, y = y1,
        ylim = c(0,max(range1,range2,range3)),
        col = 'red',
        lwd = line_width,
        type = 'l',
        main = "probability density of three clusters",
        xlab = "x",
        ylab = "density")
    }
    lines(x = x, y = y2, col = 'blue', lwd = line_width)
    lines(x = x, y = y3, col = 'black', lwd = line_width)
    lines(x = x, y = y1+y2+y3, col = 'brown', lwd = line_width)
  }

  legend(x="topleft",legend=c("cluster 1", "cluster 2", "cluster 3", "summation"), col=c("red", "blue", "black", "brown"), lty=rep(1,4), lwd = line_width)

}
