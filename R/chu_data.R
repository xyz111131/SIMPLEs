#' Single cell RNASeq and bulk RNASeq data of human embryonic stem cell differentiation.
#'
#' Data from a study of human embryonic stem cell differentiation towards definitive endoderm. It includes both the single cell RNASeq and the bulk RNASeq data. Here we randomly choose part of the data.
#'
#' @docType data
#'
#' @usage data(chu)
#'
#' @format include chu_normalized_data, chu_bulk_mean, differential_genes and chu_cell_type randomly selected from the raw data.
#'
#' @keywords datasets
#'
#' @references Chu, L.F. et al., (2016) Genome Biology 17:1-20
#' (\href{https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1033-x}{Genome Biology})
#'
#' @source \href{https://phenome.jax.org/projects/Moore1b}{QTL Archive}
#'
#' @examples
#' data(chu)
#' genes <- nrow(chu_normalized_data)
#' cells <- ncol(chu_normalized_data)
#' stat_celltype <- table(factor(chu_cell_type))
"chu"
