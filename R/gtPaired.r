# created on Jan 6, 2017
# globaltest wrapper

gtPaired.default<-function(y, alternative, null, phenoDat,
  verbose = FALSE)
{

  ans<-try(globaltest::gt(response=y, alternative=alternative, 
    null=null, data=phenoDat),
    silent=TRUE)

  aaa<-attr(ans,which="class")
  if(length(aaa)>1)
  { aaa<-aaa[1] }
  if (aaa == "try-error")
  {
    pval <- NA
    stat <- NA
    cat("try-error; set pval=NA and stat=NA\n")
  } else {

    if(verbose)
    { print(ans) }

    stat = z.score(ans)
    pval = p.value(ans)
  }

  res<-c(stat, pval)
  names(res)<-c("stat", "pval")
  return(res)

}

gtPaired<-function(es, 
  alpha=0.05, pvalAdjMethod="fdr")
{
  dat=exprs(es)
  phenoDat=pData(es)


  nPairs = ncol(dat)
  alternative<-matrix(rep(1, nPairs), ncol=1)
  colnames(alternative)<-"intercept"
  null = ~0

  mat<-t(apply(dat, 1, gtPaired.default,
    alternative=alternative, null=null, 
    phenoDat = phenoDat, verbose=FALSE))
  colnames(mat)<-c("stat", "pval")

  stat<-as.numeric(mat[,1])
  pval<-as.numeric(mat[,2])

  # significant results are always positive based on
  # the special setting of alternative and null
  # so we manually set sign for stat based on mean of log2 difference
  m=apply(dat, 1, mean, na.rm=TRUE)
  neg=which(m<0)
  stat[neg] = -stat[neg]

  # get significant gene list with FDR adjusted p.values less than 0.01
  p.adj<-p.adjust(pval, method=pvalAdjMethod)

  resFrame=data.frame(probeid=featureNames(es), stat=stat, pval=pval,
    p.adj=p.adj)

  # the "higher" value of memSubj is over-expressed
  pos1<-which(p.adj<alpha & stat>0)
  len1<-length(pos1)

  # the "higher" value of memSubj is under-expressed
  pos2<-which(p.adj<alpha & stat<0)
  len2<-length(pos2)

  nGenes<-nrow(dat)
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

