\name{gtPaired}
\alias{gtPaired}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Wrapper function for gt
}
\description{
Wrapper function for the function gt in the Bioconductor package globaltest.
}
\usage{
gtPaired(
  es, 
  alpha = 0.05, 
  pvalAdjMethod = "fdr")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{es}{
An \code{ExpressionSet} object stores within-pair log2 difference.
}
  \item{alpha}{
cutoff for adjusted p-value. If an adjusted pvalue is less than
\code{alpha}, we claim this probe is significant.
}
  \item{pvalAdjMethod}{
character. indicating which p-value adjustment method will be used.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
A list with 3 elements:
\item{resFrame }{a data frame with 4 columns: probeid, stat, pval, and p.adj.}
\item{memGenes }{3-cluster probe cluster membership. 1 indicates over-expressed probes;
 2 indicates under-expressed probes; 3 indicates non-differentially expressed
probes.}
\item{memGenes2 }{2-cluster probe cluster membership. 1 indicates differentially-expressed probes;
 0 indicates non-differentially expressed probes.}

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

result = gtPaired(es = E_Set)

print(table(result$memGenes, fData(E_Set)$memGenes.true))
print(table(result$memGenes2, fData(E_Set)$memGenes2.true))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }% use one of  RShowDoc("KEYWORDS")

