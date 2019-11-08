\name{scimpclu}
\alias{scimpclu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Impute zero entries and clustering for scRNASeq data
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Impute zeros in the gene expression data using the expression level in similar cells and gene-gene correlation. Zero entries in the observed expression matrix come from molecule loss during the experiment ('dropout') or too low expression to be measured. We used Monte Carlo EM algorithm to sample the imputed values and learn the parameters.

}
\usage{
scimpclu(dat, K0, M0 = 1, iter=10, est_lam = 2, impt_it = 1, penl =1 ,sigma0 = 100,  pi_alpha = 1, beta = NULL, verbose=F, max_lambda=F, lambda =NULL, sigma=NULL, mu=NULL, est_z = 1, clus = NULL,p_min = 0.8, cutoff=0.1, K = 10, min_gene = 1000, num_mc = 3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{dat}{scRNASeq data matrix. Each row is a gene, each column is a cell.}
  \item{K0}{Number of latent gene modules. See details.}
  \item{M0}{Number of clusters. See details.}
  \item{clus}{Initial clustering of scRNASeq data. If NULL, the function will use PCA and Kmeans to do clustering initially.}
  \item{K}{The number of PCs used in the initial clustering. Default = 20. }
  \item{iter}{Number of EM iterations for full data set. See details. }
  \item{est_z}{The iteration starts to update z.}
  \item{impt_it}{The iteration starts to sample new imputed values in initial phase. See details.}
   \item{max_lambda}{Whether to maximize over lambda.}
   \item{est_lam}{The iteration starts to estimate lambda.}
  \item{penl}{L1 penalty for the factor loadings.}
  \item{sigma0}{The variance of the prior distribution of \eqn{\mu}.}
  \item{pi_alpha}{The hyperparameter of the prior distribution of \eqn{\pi}. See details. }
   \item{beta}{A G by K0 matrix. Initial values for factor loadings (B). If null, beta will initialze from normal distribution with mean zero and variance M0/K0. See details. }
    \item{lambda}{A M0 by K0 matrix. Initial values for the variances of factors. Each column is for a cell cluster. If null, lambda will initialize to be 1/M0. See details. }
     \item{sigma}{A G by M0 matrix. Initial values for the variance of idiosyncratic noises. Each column is for a cell cluster. If null, sigma will initialize to be 1. See details. }
     \item{mu}{A G by M0 matrix. Initial values for the gene expression mean of each cluster. Each column is for a cell cluster. If NULL, it will take the sample mean of cells weighted by the probability in each cluster. See details.}
    \item{p_min}{Initialize parameters using genes expressed in at least \emph{p_min} proportion of cells. If the number genes selected is less than \emph{min_gene}, select \emph{min_gene} genes with higest proportion of non zeros. Default = 0.8. }
    \item{min_gene}{Minimal number of genes used in the initial phase. See details. }
    \item{cutoff}{The value below cutoff is treated as no expression. Default = 0.1. }
   \item{verbose}{Whether to show some intermediate results. Default = False.}
   \item{num_mc}{The number of Gibbs steps to generate new imputed data.}
}
\details{
We assume that the cells come from M0 clusters. Within each cell cluster, the 'true' gene expression is modeled by a multivariate Gaussian distribution whose covariance matrix can be composed into a low rank matrix (a couple of latent gene modules) and idiosyncratic noises. Gene modules are shared among cell clusters though the coexpression level of each gene module can be different.

Suppose there are G genes and n cells. For each cell cluster, the gene expression follows \eqn{Y|Z=m~MVN(\mu_m, B\Lambda_m B^T + \Sigma_m)} where B is a G by K0 matrix, \eqn{\Sigma_m} is a G by G diagonal matrix whose diagonal entries are specified by \emph{sigma}, and \eqn{\Lambda_m} is a K0 by K0 diagonal matrix whose diagonal entries are specified by \emph{lambda}. \eqn{P(Z_m) = \pi_m} where \eqn{\pi~Dir(\alpha)}.

The algorithm first runs Monte Carlo EM using only the genes with low dropout rate (initial phase) and initializes factor loadings and clustering membership. Then it runs another rounds of Monte Carlo EM using all the genes. In the initial phase, we used the genes with dropout rate less than \emph{1 - p_min}; if the number of genes is less than \emph{min_gene}, we ranked the genes by the number cells with nonzero expression and kept the top \emph{min_gene} genes.
}
\value{
\item{loglik}{The log-likelihood of the imputed gene expression at each iteration.}
\item{pi}{Probabilites of cells belong to each cluster.}
\item{mu}{Mean expression for each cluster}
\item{sigma}{Variances of idiosyncratic noises for each cluster.}
\item{beta}{Factor loadings.}
\item{lambda}{Variances of factors for each cluster.}
\item{z}{The probability of each cell belonging to each cluster.}
\item{Ef}{Conditonal expection the factors for each cluster \eqn{E(f_i|z_i = m)}. A list with length M0, each element in the list is a n by K0 matrix.}
\item{Varf}{Conditonal covariance of factors for each cluster \eqn{Var(f_i|z_i = m)}. A list with length M0, each element in the list is a K0 by K0 matrix.}
\item{Yimp0}{A matrix contains the expectation of imputed expression.}
\item{Y}{Last sample of imputed matrix.}
\item{pg}{A G by M0 matrix, dropout rate for each gene in each cluster.}
\item{geneM}{Gene mean. If centerized each gene before estimating the parameters, provide the overall mean of gene expression removed from the data matrix. }
\item{geneSd}{Gene standard deviation. If scaled each gene before estimating the parameters, provide the overall standard deviation of gene expression removed from the data matrix. }
\item{initclus}{Output initial cluster results.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Z. Hu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
scimpclu_bulk
}
\examples{

library(foreach)
library(doParallel)
library(SIMPLE)

source("SIMPLE/utils.R")
M0 = 3  # simulate number of clusters
n = 300 # number of cells

simu_data = simulation_bulk(n=300, S0 = 20, K = 6, MC=M0, block_size = 32, indepG = 1000 - 32*6, verbose=F, overlap=0)

Y2 = simu_data$Y2
K0 = 6 # number of factors

registerDoParallel(cores = 6)  # parallel

# estimate the parameters
result <- scimpclu(Y2, K0, M0, celltype=rep(1, n), clus = NULL, K = 20, p_min = 0.5, max_lambda=T, min_gene = 200,cutoff=0.01)

# sample imputed values
result2 = do_impute(Y2, result$Y, result$beta, result$lambda, result$sigma, result$mu, result$pi, result$geneM, result$geneSd, rep(1, n), mcmc=50, burnin = 5, pg = result$pg, cutoff = 0.01)

# evaluate cluster performance
celltype_true = simu_data$Z
mclust::adjustedRandIndex(apply(result$z,1, which.max), celltype_true)
# or redo clustering based on imputed values (sometimes work better for real data)
getCluster(result2$impt, celltype_true, Ks = 20, M0 = M0)[[1]]

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
