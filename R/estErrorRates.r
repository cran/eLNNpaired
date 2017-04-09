# created on Jan. 4, 2017
#  estimate error rates: FDR, FNDR, FPR, FNR

# memGenes2.est - estimated cluster membership
#                 memGenes2.est=1 means differentially expressed
#                 memGenes2.est=0 means non-differentially expressed
# memGenes2.est - estimated cluster membership
#                 memGenes2.est=1 means differentially expressed
#                 memGenes2.est=0 means non-differentially expressed
estErrorRates = function(memGenes2.est, memGenes2.true)
{
  # FDR
  pos.numer = which(memGenes2.est==1 & memGenes2.true==0)
  pos.denom = which(memGenes2.est==1)

  if(length(pos.denom))
  {
    FDR = length(pos.numer)/length(pos.denom)
  } else {
    FDR = 0
  }

  # FNDR
  pos.numer = which(memGenes2.est==0 & memGenes2.true==1)
  pos.denom = which(memGenes2.est==0)

  if(length(pos.denom))
  {
    FNDR = length(pos.numer)/length(pos.denom)
  } else {
    FNDR = 0
  }

  # FPR
  pos.numer = which(memGenes2.est==1 & memGenes2.true==0)
  pos.denom = which(memGenes2.true==0)

  if(length(pos.denom))
  {
    FPR = length(pos.numer)/length(pos.denom)
  } else {
    FPR = 0
  }

  # FNR
  pos.numer = which(memGenes2.est==0 & memGenes2.true==1)
  pos.denom = which(memGenes2.true==1)

  if(length(pos.denom))
  {
    FNR = length(pos.numer)/length(pos.denom)
  } else {
    FNR = 0
  }

  res=c(FDR, FNDR, FPR, FNR)
  names(res) = c("FDR", "FNDR", "FPR", "FNR")

  return(res)
}

