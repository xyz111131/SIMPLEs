\name{do_impute}
\alias{do_impute}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sampling the imputed values and factors at given parameters.
}
\description{
Get the posterior mean and variance for the imputed values and factors given the parameters, which can be returned from \emph{scimpclu_bulk}.
}
\usage{
do_impute(dat, Y, beta, lambda, sigma, mu, pi, pos_mean = NULL, pos_sd = NULL, celltype = NULL, mcmc=10, burnin = 2, verbose=F, pg = 0.5,cutoff=0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{scRNASeq data matrix. Each row is a gene, each column is a cell.}
  \item{Y}{Initial imputed data. Can return from \emph{scimpclu_bulk}.}
  \item{beta}{Factor loadings.}
  \item{lambda}{Variances of factors for each cluster.}
  \item{sigma}{Variances of idiosyncratic noises for each cluster.}
  \item{mu}{Mean expression for each cluster}
  \item{pi}{Probabilites of cells belong to each cluster.}
  \item{pos_mean}{Gene mean. If centerized each gene before estimating the parameters, provide the overall mean of gene expression removed from the data matrix. }
  \item{pos_sd}{Gene standard deviation. If scaled each gene before estimating the parameters, provide the overall standard deviation of gene expression removed from the data matrix.}

}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{loglik}{The log-likelihood of the imputed gene expression at each iteration.}
\item{impt}{A matrix contains the expectation of imputed expression.}
\item{EF}{Posterior means of factors}
\item{varF0}{Posterior mean of the covariance between factors}
\item{varF}{Posterior variances of factors}
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
scimpclu_bulk
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x)
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
