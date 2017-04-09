
#' Generate an ExpressionSet from eLNN model.
#'
#' @return An ExpressionSet with exprs() from our model and fData$true_cluster filled by genes' true cluster information.
#' @param G An integer, the number of genes.
#' @param n An integer, the number of pairs for each gene.
#' @param psi A vector of length 10. It contains the parameters after reparameterization as illustrated in paper.
#' @param t_pi A vector of length 3 and summation equal to 1. It is 'pi' as in paper, we do not use 'pi' since it is reserved as a mathematical constatnt in R.
#' @param c1 A parameter in constraints. It should be in the form of c1 = qnorm(X), where X is a decimal smaller than 1 but close to 1. Larger X gives more stringent constraint. Default value is c1 = qnorm(0.95).
#' @param c2 A parameter in constraints. It should be in the form of c2 = qnorm(Y), where Y is a decimal larger than 0 but close to 0. Smaller Y gives more stringent constraint. Default value is c2 = qnorm(0.05).
#' @export
#' @example example/example.r
#gen_eLNNpaired <- function(
#  G,
#  n,
#  psi,
#  t_pi,
#  c1 = qnorm(0.95),
#  c2 = qnorm(0.05))
#{
#  dat = matrix(, nrow = G, ncol = n)
#  cluster_info = matrix(rep(0,G*3),G,3)
#
#  colnames(cluster_info) = c("true_cluster","est_cluster","flag")
#
#  for (row in 1:G)
#  {
#    temp = runif(1)
#    if (temp<t_pi[1]) cluster = 1
#    else if(temp<t_pi[1]+t_pi[2]) cluster = 2
#    else cluster = 3
#
#    cluster_info[row,1] = cluster
#
#    if (cluster == 1)
#      {
#        mu_0 = exp(psi[cluster*4-3])
#        k = pnorm(psi[cluster*4-2])
#        beta = exp(psi[cluster*4])
#        alpha = exp(psi[cluster*4-1]) + 1 + beta * c1^2 * (1+sqrt(k))^2/mu_0^2
#      }
#      else if (cluster == 2)
#      {
#        mu_0 = -exp(psi[cluster*4-3])
#        k = pnorm(psi[cluster*4-2])
#        beta = exp(psi[cluster*4])
#        alpha = exp(psi[cluster*4-1]) + 1 + beta * c2^2 * (1+sqrt(k))^2/mu_0^2
#      }
#      else
#      {
#        mu_0 = 0
#        k = pnorm(psi[cluster*4-2])
#        beta = exp(psi[cluster*4])
#        alpha = exp(psi[cluster*4-1])
#      }
#
#    tau_g = rgamma(1, alpha, beta)
#
#    if (cluster == 3)
#    {
#      mu_g = 0
#    }
#    else{
#      mu_g = rnorm(1, mean = mu_0, sd = sqrt(k/tau_g))
#    }
#
#    dat[row,] = rnorm(n,mean = mu_g, sd = sqrt(1/tau_g))
#  }
#
#  return (ExpressionSet(assayData = dat, featureData = new("AnnotatedDataFrame", data = data.frame(cluster_info))))
#}
#

# psi is 10x1 vector of reparameterized parameters
#  c(delta1, xi1, lambda1, nu1,
#    delta2, xi2, lambda2, nu2,
#    lambda3, nu3)
#
# t_pi = c(pi_1, pi_2)
gen_eLNNpaired <- function(
  G,
  n,
  psi,
  t_pi,
  c1 = qnorm(0.95),
  c2 = qnorm(0.05))
{
  dat = matrix(NA, nrow = G, ncol = n)
  pi1=t_pi[1]
  pi2=t_pi[2]
  pi3=1-pi1-pi2

  memGenes=sample(x=c(1,2,3), size=G, prob=c(pi1, pi2, pi3), replace=TRUE)
  memGenes2=rep(1, G)
  pos=which(memGenes==3)
  if(length(pos))
  {	   
    memGenes2[pos]=0
  }

  pos1=which(memGenes==1)
  pos2=which(memGenes==2)
  pos3=which(memGenes==3)

  if(length(pos1) == 0 || length(pos2) == 0 || length(pos3)==0)
  {
    stop("each probe cluster should have at least 2 probes!\nPlease re-assign values of t_pi!\n")
  }

  G1=length(pos1)
  G2=length(pos2)
  G3=length(pos3)

  # get original-scale parameters
  para.orig=getPara.orig(
    delta1=psi[1], xi1=psi[2], lambda1=psi[3], nu1=psi[4],
    delta2=psi[5], xi2=psi[6], lambda2=psi[7], nu2=psi[8],
    lambda3=psi[9], nu3=psi[10],
    c1=c1, 
    c2=c2
  )

  # cluster 1
  mu_0 = para.orig[1]
  k = para.orig[2]
  alpha = para.orig[3]
  beta = para.orig[4]
  tau_g = rgamma(G1, alpha, beta)
  for(g in 1:G1)
  {
    mu_g = rnorm(1, mean = mu_0, sd = sqrt(k/tau_g[g]))
    dat[pos1[g],] = rnorm(n, mean = mu_g, sd = sqrt(1/tau_g[g]))
  }

  # cluster 2
  mu_0 = para.orig[5]
  k = para.orig[6]
  alpha = para.orig[7]
  beta = para.orig[8]
  tau_g = rgamma(G2, alpha, beta)
  for(g in 1:G2)
  {
    mu_g = rnorm(1, mean = mu_0, sd = sqrt(k/tau_g[g]))
    dat[pos2[g],] = rnorm(n, mean = mu_g, sd = sqrt(1/tau_g[g]))
  }

  # cluster 3
  alpha = para.orig[9]
  beta = para.orig[10]
  tau_g = rgamma(G3, alpha, beta)
  for(g in 1:G3)
  {
    dat[pos3[g],] = rnorm(n, mean = 0, sd = sqrt(1/tau_g[g]))
  }

  probeid=1:G
  subjid=paste("subj", 1:n, sep="")

  rownames(dat) =probeid
  colnames(dat)= subjid

  pDat=data.frame(subjid=subjid)
  rownames(pDat)=subjid

  gene=paste("g", 1:G, sep="")
  fDat=data.frame(probeid=probeid, gene=gene, chr=rep(1, G),
    memGenes.true=memGenes, memGenes2.true=memGenes2)
  rownames(fDat)=probeid

  es=iCheck::genExprSet(
       ex=dat,
       pDat=pDat,
       fDat = fDat,
       annotation = "")

  invisible(es)
}

