\name{estErrorRates}
\alias{estErrorRates}
\title{
Estimate Error Rates
}
\description{
Estimate error rates: false discovery rate (FDR), 
false non-discovery rate (FNDR), false positve rate (FPR), 
and false negative rate (FNR) for data with known
cluster membership for gene probes.
}
\usage{
estErrorRates(memGenes2.est, memGenes2.true)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{memGenes2.est}{
estimated cluster membership for gene probes. \code{memGenes2.est}
must take only two values 1 and 0.
}
  \item{memGenes2.true}{
true cluster membership for gene probes. \code{memGenes2.true}
must take only two values 1 and 0.
}
}
\details{
FDR is the proportion of detected differentially expressed (DE) gene probes
that are truly non-differentially expressed (NE) gene probes
among all detected DE gene probes.

FNDR is the proportion of detected NE gene probes
that are truly DE gene probes
among all detected NE gene probes.

FPR is the proportion of detected DE gene probes
that are truly NE gene probes
among all truly NE gene probes.

FNR is the proportion of detected NE gene probes
that are truly DE gene probes
among all truly DE gene probes.

}
\value{
A vector of 4 elements: FDR, FNDR, FPR, and FNR.
}
\references{Li Y, Morrow J, Raby B, Tantisira K, Weiss ST, Huang W, Qiu W. (2017), <doi:10.1371/journal.pone.0174602>}
\author{
Yunfeng Li <colinlee1999@gmail.com> and Weiliang Qiu <stwxq@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
  memGenes2.est = c(1, 1, 1, 0, 0, 0)
  memGenes2.true = c(1, 0, 1, 0, 0, 0)

  err = estErrorRates(memGenes2.est = memGenes2.est, 
    memGenes2.true = memGenes2.true)
  print(err)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }% use one of  RShowDoc("KEYWORDS")

