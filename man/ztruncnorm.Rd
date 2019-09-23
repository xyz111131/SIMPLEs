\name{ztruncnorm}
\alias{ztruncnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a zero-inflated censored normal distribution.
}
\description{
Fit each gene expression by a zero-inflated censored normal distribution by direct maximize the log-likelihood and return MLE of mean, variance and probability of mixture component.
}
\usage{
ztruncnorm(dat, cutoff=0.5, iter = 100, s_upper = 4, s_lower = 0.1, p_min = 0.6, p_max=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{Data matrix.}
  \item{cutoff}{The value below cutoff is treated as no expression. Default = 0.1. }
  \item{iter}{Max number of iteration for optimization. Default = 100.}
  \item{s_upper}{upper bound of standard deviation. Default = 4. }
  \item{s_lower}{lower bound of standard deviation. Default = 0.1.}
  \item{p_min}{lower bound of the probability of the normal component. Default = 0.6. }
  \item{p_max}{upper bound of the probability of the normal component. Can be either a scaler or a vector. If it is a vector, each element is the upper bound for a gene. Default = 1. }
}
\details{
It will return mean -1, std 0.5 and p = p_max for the genes that the number of nonzero entries is less than 2 to indicate that the gene has no expression and will return NA if the optimization is not convergence or error which will be filled with sample mean and std, p_max afterwards.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{param}{A matrix with mean and standard deviation for each gene. }
\item{pg}{A vector for the probablity of the normal component.}
%% ...
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
}

