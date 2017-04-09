\name{samrOneClass}
\alias{samrOneClass}
\title{
A Wrapper Function of samr for Paired Data
}
\description{
A wrapper function of samr for paired data.
}
\usage{
samrOneClass(
  es, 
  fdr.output = 0.05, 
  nperms = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{es}{
An \code{ExpressionSet} object stores within-pair log2 difference.
}
  \item{fdr.output}{
fdr cutoff for output in
          significant genes table. (see \code{SAM} function in
\code{samr} package.
}
  \item{nperms}{
Number of permutations used to estimate false discovery rates.
(see \code{samr} function in \code{samr} package).
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list of 5 elements:
\item{samr.obj }{object returned by \code{samr} function}
\item{del}{estimated cutoff value required by \code{samr.compute.siggenes.table} function. See the source code of \code{SAM} function}
\item{delta.table}{object returned by \code{samr.compute.delta.table}}
\item{siggenes.table}{object returned by \code{samr.compute.siggenes.table} function}
\item{memGenes}{probe cluster membership. 1 indicates over-expressed probes;
2 indicates under-expressed probes; 3 indicates non-differentially
expressed probes}
\item{memGenes2}{2-cluster probe cluster membership. 1 indicates 
differentially-expressed probes;
0 indicates non-differentially-expressed probes.}

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
set.seed(100)
G = 500
n = 10

delta_1 = -0.8184384  
xi_1 = -1.1858546 
lambda_1 = -10.6309216  
nu_1 = -3.5536255  

delta_2 = -0.8153614  
xi_2 = -1.4120148 
lambda_2 = -13.1999427  
nu_2 = -3.3873531   

lambda_3 = 0.7597441  
nu_3 = -2.0361091 

psi = c(delta_1, xi_1, lambda_1, nu_1,
        delta_2, xi_2, lambda_2, nu_2,
        lambda_3, nu_3)
t_pi = c(0.08592752, 0.07110449)

c1 = qnorm(0.95)
c2 = qnorm(0.05)

E_Set = gen_eLNNpaired(G, n, psi, t_pi, c1, c2)

result = samrOneClass(es=E_Set)

print(table(result$memGenes, fData(E_Set)$memGenes.true))
print(table(result$memGenes2, fData(E_Set)$memGenes2.true))

}
\keyword{ method }% use one of  RShowDoc("KEYWORDS")

