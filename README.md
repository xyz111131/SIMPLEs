# SIMPLE: SIngle cell RNASeq iMPutation and celL clustEring
=====================================
Zhirui Hu

## Introduction
SIMPLE is a R package that imputes "dropout" values (i.e. zero or small entries)
in the single cell RNASeq data based on cell similarities and gene correlations.
The imputed matrix can be used for dimension reduction and visualization of the
cell spectrum, to identify markers between different samples and construct gene
co-expression network. SIMPLEs iteratively clusters the cells, identifies
correlated gene modules infers the probability of the dropout event for each
zero entry, imputes zeros within each cluster utilizing the expression of other
correlated genes. It will impute the technical zeros (or "dropout") while
maintain the biological zeros at low expression level.

The imputation process is based on the correlations between genes within similar
cell types, which is modeled by several common latent factors or gene modules as
well as the gene-specific dropout rate. Although the dropout rate can be
estimated from the empirical distribution of gene expression in the scRNASeq, it
could interference with estimating the gene correlation structure, especially
for lowly expressed genes. Integrating with the corresponding bulk RNASeq data
which serves as the average gene expression across cells, provides extra sources
of information on the dropout rate per gene. It can give a better estimate of
the "dropout" rate which would influence how much the data should be imputed. We
called our method integrating bulk RNASeq data as **SIMPLE-B**, otherwise
**SIMPLE**. We referred our toolset including SIMPLE and SIMPLE-B as
**SIMPLEs**.

## Installation
------------
The package is not on CRAN yet. You can use the following codes in `R` to
install it.

```r
install.packages("devtools")
library(devtools)_
devtools::install_github("xyz111131/SIMPLEs")
```
