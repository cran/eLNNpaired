# created on Jan. 3, 2017 by Weiliang Qiu

# wrapper for one sample test using samr
samrOneClass = function(es, fdr.output=0.05, nperms=100)
{

  x=exprs(es)
  probes=1:nrow(x)
  genes=1:nrow(x)
  mydat=list(x=x,y=rep(1, ncol(x)), geneid=probes,
     genenames=genes, logged2=TRUE)

  # prevent output "perm=.." to screen
  mylog <- capture.output({
    samr.obj<-samr(data=mydat,  resp.type="One class", nperms=nperms)
  })


  # prevent output "perm=.." to screen
  mylog2 <- capture.output({
    delta.table <- samr.compute.delta.table(samr.obj)
  })

  siggenes.table <- del <- NULL
  delta.table <- delta.table[delta.table[, "# called"] > 0,
      , drop = FALSE]
  memGenes=rep(3, nrow(x))

  if (nrow(delta.table) > 0) {
      oo <- which(delta.table[, "median FDR"] >= fdr.output)
      if (length(oo) > 0) {
          oo <- oo[length(oo)]
      }
      else {
          oo <- 1
      }
      delta.table <- delta.table[oo:nrow(delta.table), , drop = FALSE]
      del <- delta.table[1, "delta"]

      # prevent output "perm=.." to screen
      mylog3 <- capture.output({
        siggenes.table <- samr.compute.siggenes.table(samr.obj=samr.obj,
          del=del, data=mydat, delta.table=delta.table)
      })

      if(siggenes.table$ngenes.up)
      {
        pos.up=as.numeric(siggenes.table$genes.up[,2])
        memGenes[pos.up]=1
      }
      if(siggenes.table$ngenes.lo)
      {
        pos.lo=as.numeric(siggenes.table$genes.lo[,2])
        memGenes[pos.lo]=2
      }
       
      siggenes.table$genes.up = siggenes.table$genes.up[, -1]
      siggenes.table$genes.lo = siggenes.table$genes.lo[, -1]

      rang = 3:6
      siggenes.table$genes.up[, rang] = as.numeric(siggenes.table$genes.up[, rang])
      siggenes.table$genes.lo[, rang] = as.numeric(siggenes.table$genes.lo[, rang])
  }
  
  memGenes2=rep(1, nrow(x))
  pos=which(memGenes==3)
  if(length(pos))
  {
    memGenes2[pos]=0
  }

  out = list(samr.obj = samr.obj, del = del, delta.table = delta.table,
      siggenes.table = siggenes.table, memGenes=memGenes,
      memGenes2=memGenes2)

  invisible(out)
}

