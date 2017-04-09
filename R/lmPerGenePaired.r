# created on Jan 6, 2017

# get two-sided pvalue based on z-value
calPval<-function(Z)
{
  # pvalue
  if(Z<=0)
  {
    pval<-2*pnorm(Z)
  } else {
    pval<-2*(1-pnorm(Z))
  }
  pval<-min(pval, 0.999999)
  return(pval)
}


# gene-wise linear model for log2 differences
lmPerGenePaired<-function(
  es,
  alpha=0.05, 
  pvalAdjMethod = "fdr"
)
{
  fit = GSEAlm::lmPerGene(eSet=es, formula=~1)

  stat<-as.vector(fit$tstat[1,])
  len<-length(stat)
  pval<-rep(0, len)
  for(i in 1:len)
  { pval[i] <- calPval(stat[i]) }

  # get significant gene list with FDR adjusted p.values less than alpha
  p.adj<-p.adjust(pval, method=pvalAdjMethod)
  resFrame=data.frame(probeid=featureNames(es), stat=stat, pval=pval,
    p.adj=p.adj)

  # the "higher" value of memSubj is over-expressed
  pos1<-which(p.adj<alpha & stat>0)
  len1<-length(pos1)

  # the "higher" value of memSubj is under-expressed
  pos2<-which(p.adj<alpha & stat<0)
  len2<-length(pos2)

  nGenes<-nrow(es)
  memGenes<-rep(3, nGenes)
  memGenes2<-rep(1, nGenes)

  # mixing proportion
  if(len1)
  {
    memGenes[pos1]<-1
  } 
  if(len2)
  {
    memGenes[pos2]<-2
  } 
  pos3=which(memGenes==3)
  if(length(pos3))
  {
    memGenes2[pos3]=0
  }

  res<-list(resFrame=resFrame, memGenes=memGenes, memGenes2=memGenes2)
  invisible(res)
}


