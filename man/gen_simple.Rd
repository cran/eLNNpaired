\name{gen_simple}
\alias{gen_simple}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate Simulated Data Sets from a Simple Model
}
\description{
Generate Simulated Data Sets from a Simple Model.
}
\usage{
gen_simple(
  G, 
  n = 30, 
  psi = c(0.441, 1, -0.442, 1, 2), 
  t_pi = c(0.086, 0.071)
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{G}{An integer, the number of genes.}

\item{n}{An integer, the number of pairs for each gene.}

  \item{psi}{
A vector of 5 elements containing model parameters
\eqn{\mu_1}{mu_1}, \eqn{\sigma_1}{sigma_1},
\eqn{\mu_2}{mu_2}, \eqn{\sigma_2}{sigma_2},
and \eqn{\sigma_3}{sigma_3}.
}
\item{t_pi}{the cluster proportion for cluster 1 (over-expressed probes)
and cluster 2 (under-expressed probes).}
}
\details{
We assume there are three clusters of gene probes: (1) over-expressed;
(2) under-expressed; and (3) non-differentially expressed.
For probes in cluster 1, we assume the within-pair log2 difference of 
gene expression is from \eqn{N(\mu_1, \sigma_1^2)}{N(mu_1, sigma_1^2)}.
For probes in cluster 2, we assume the within-pair log2 difference of 
gene expression is from \eqn{N(\mu_2, \sigma_2^2)}{N(mu_2, sigma_2^2)}.
For probes in cluster 3, we assume the within-pair log2 difference of 
gene expression is from \eqn{N(0, \sigma_3^2)}{N(0, sigma_3^2)}.
\eqn{\mu_1>0}{mu_1>0} and \eqn{\mu_2<0}{mu_2<0}.
}
\value{
An ExpressionSet object, the feature data frame of which include
\code{memGenes.true} (3-cluster membership for gene probes)
and \code{memGenes2.true} (2-cluster membership for gene probes).

In 3-cluster membership, 1 indicates over-expressed,
2 indicates under-expressed, and 3 indicates non-differentially expressed.

In 2-cluster membership, 1 indicates differentially expressed,
0 indicates non-differentially expressed.

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

es=gen_simple(
  G = 500,
  n = 30,
  psi = c(0.441, 1, -0.442, 1, 2),
  t_pi = c(0.086, 0.071)
)

print(es)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }% use one of  RShowDoc("KEYWORDS")

