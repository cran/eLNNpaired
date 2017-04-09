# created on Jan 6, 2017 based on file 'generate_data.r'

# generate log2 difference for each pair of samples
#
# G - number of probes
# n - number of pairs
# psi=c(mu_1, sigma_1, mu_2, sigma_2, sigma_3)
# t_pi=c(pi_1, pi_2)
#   pi_1 - proportion of over-expressed probes
#   pi_2 - proportion of under-expressed probes

gen_simple <- function(
  G,
  n = 30,
  psi = c(0.441, 1, -0.442, 1, 2), 
  t_pi = c(0.086, 0.071)
)
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

  mu_1 = psi[1] 
  sigma_1 = psi[2] 
  mu_2 = psi[3] 
  sigma_2 = psi[4] 
  mu_3 = 0 
  sigma_3 = psi[5] 
  
  over_expressed_genes = rnorm(
  	n * G1, 
  	mean = mu_1, 
  	sd = sigma_1)

  dat[pos1,] = matrix(
  	over_expressed_genes, 
  	nrow=G1, 
  	ncol=n)
  
  under_expressed_genes = rnorm(
  	n * G2,
  	mean = mu_2,
  	sd = sigma_2)

  dat[pos2,] = matrix(
  	under_expressed_genes, 
  	nrow=G2, 
  	ncol=n)
  
  non_diff_expressed_genes = rnorm(
  	n * G3,
  	mean = mu_3,
  	sd = sigma_3)

  dat[pos3,] = matrix(
  	non_diff_expressed_genes,
  	nrow=G3,
  	ncol=n)
  
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

