\name{init_impute_bulk}
\alias{init_impute_bulk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Initialize imputed matrix integrating bulk RNASeq
}
\description{

}
\usage{
init_impute_bulk(Y2, clus, bulk, pg1, cutoff = 0.1,verbose=F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y2}{ scRNASeq data matrix. Each row is a gene, each column is a cell.}
  \item{clus}{A numeric vector for the cell type labels of cells in the scRNASeq. Each cell type has corresponding the mean expression in the bulk RNASeq data. The labels must start from 1 to the number of types. See details. }
  \item{bulk}{Bulk RNASeq data matrix. Each row is a gene which must be ordered the same as the scRNASeq data. Each column is a cell type which must be ordered as the cell type label in \emph{clus}. }
  \item{pg1}{A data matrix for the dropout rate for each gene which is the ratio between mean expression in the scRNASeq and bulk RNASeq. Each row is a gene which must be ordered the same as the scRNASeq data. Each column is a cell type which must be ordered as the cell type label in \emph{clus}.}
  \item{cutoff}{The value below cutoff is treated as no expression. Default = 0.1.}
  \item{verbose}{Whether to plot some intermediate result. Default = False.}
}
\details{
If the scRNASeq data does not have the cell type information, then all cells can be treated as a single type and the input bulk RNASeq should also have one column that is the mean expression over all cell types.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
Imputed gene expression matrix
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
