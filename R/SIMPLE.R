#' Initialize imputed matrix integrating bulk RNASeq
#'
#' @details
#'If the scRNASeq data does not have the cell type information, then all cells
#'   can be treated as a single type and the input bulk RNASeq should also have
#'   one column that is the mean expression over all cell types.
#' @param Y2  scRNASeq data matrix. Each row is a gene, each column is a cell.
#' @param clus A numeric vector for the cell type labels of cells in the scRNASeq. Each cell type has corresponding the mean expression in the bulk RNASeq data. The labels must start from 1 to the number of types. See details.
#' @param bulk Bulk RNASeq data matrix. Each row is a gene which must be ordered the same as the scRNASeq data. Each column is a cell type which must be ordered as the cell type label in \emph{clus}.
#' @param pg1 A data matrix for the dropout rate for each gene which is the ratio between mean expression in the scRNASeq and bulk RNASeq. Each row is a gene which must be ordered the same as the scRNASeq data. Each column is a cell type which must be ordered as the cell type label in \emph{clus}.
#' @param cutoff The value below cutoff is treated as no expression. Default = 0.1.
#' @param verbose Whether to plot some intermediate result. Default = False.
#'
#' @return Imputed gene expression matrix
#' @author Zhirui Hu, \email{zhiruihu@g.harvard.edu}
#' @author Songpeng Zu, \email{songpengzu@g.harvard.edu}

init_impute <- function(Y2, M0, clus, p_min = 0.6, cutoff = 0.5, verbose = F) {
  # fit truncnorm for each cluster
  impute <- Y2
  G <- nrow(Y2)
  pg <- matrix(0, G, M0)
  for (i in 1:M0)
  {
    temp_dat <- as.matrix(Y2[, clus == i])
    result <- ztruncnorm(temp_dat, cutoff = cutoff, p_min = p_min) # , s_upper=3
    mu <- result[[1]][, 1]
    sd <- result[[1]][, 2]
    pg[, i] <- result[[2]]

    pg[pg > 0.99] <- 0.99

    ind <- which(is.na(mu))
    if (verbose) print(paste(length(ind), "not fit TN"))
    mu[ind] <- rowMeans(Y2[ind, clus == i, drop = F])
    sd[ind] <- rowSds(Y2[ind, clus == i, drop = F])
    pg[ind, i] <- 1

    if (verbose) print(summary(sd))

    for (j in which(clus == i))
    {
      ind <- which(Y2[, j] <= cutoff)
      # impute Y2
      ms <- mu[ind]
      sds <- sd[ind]
      p <- pg[ind, i]
      prob <- pnorm(cutoff, mean = ms, sd = sds) # compute x<0 prob
      prob_drop <- (1 - p) / (prob * p + (1 - p))
      I_drop <- rbinom(length(prob_drop), 1, prob_drop)

      # imputation for dropout
      impute[ind[I_drop == 1], j] <- rnorm(sum(I_drop == 1), ms[I_drop == 1], sds[I_drop == 1])
      # åimpute[ind[I_drop==1],j] = ms[I_drop==1]

      # imputation for non-dropout
      if (sum(I_drop == 0) > 0) {
        # r = (cutoff - ms[I_drop==0])/sds[I_drop==0]
        impute[ind[I_drop == 0], j] <- rtnorm(sum(I_drop == 0), upper = cutoff, mean = ms[I_drop == 0], sd = sds[I_drop == 0]) # ms[I_drop==0] - sds[I_drop==0] * dnorm(r)/pnorm(r)
      }
    }
  }

  impute[is.na(impute)] <- 0

  # print(sum(impute < -2))
  if (verbose) print(summary(c(impute[Y2 == 0])))
  # impute[impute < -3] = -3

  return(list(impute, pg))
}

#' SIMPLE: Imputing zero entries and clustering for scRNASeq data.
#'
#' \code{SIMPLE} imputes zeros in the gene expression data using the expression level in
#' similar cells and gene-gene correlation. Zero entries in the observed expression matrix
#' come from molecule loss during the experiment ('dropout') or too low expression to
#' be measured. We used Monte Carlo EM algorithm to sample the imputed values and
#' learn the parameters.
#'
#' @details
#' We assume that the cells come from M0 clusters. Within each cell
#' cluster, the 'true' gene expression is modeled by a multivariate Gaussian
#' distribution whose covariance matrix can be composed into a low rank matrix
#' (a couple of latent gene modules) and idiosyncratic noises. Gene modules are
#' shared among cell clusters though the coexpression level of each gene module
#' can be different. \cr
#' Suppose there are G genes and n cells. For each cell
#' cluster, the gene expression follows \eqn{Y|Z=m~MVN(\mu_m, B\Lambda_m B^T +
#' \Sigma_m)} where B is a G by K0 matrix, \eqn{\Sigma_m} is a G by G diagonal
#' matrix whose diagonal entries are specified by \emph{sigma}, and
#' \eqn{\Lambda_m} is a K0 by K0 diagonal matrix whose diagonal entries are
#' specified by \emph{lambda}. \eqn{P(Z_m) = \pi_m} where \eqn{\pi~Dir(\alpha)}. \cr
#'
#' The algorithm first runs Monte Carlo EM using only the genes with low dropout
#' rate (initial phase) and initializes factor loadings and clustering
#' membership. Then it runs another rounds of Monte Carlo EM using all the
#' genes. In the initial phase, we used the genes with dropout rate less than
#' \emph{1 - p_min}; if the number of genes is less than \emph{min_gene}, we
#' ranked the genes by the number cells with nonzero expression and kept the top
#' \emph{min_gene} genes.
#'
#' @param dat scRNASeq data matrix. Each row is a gene, each column is a cell.
#' @param K0 Number of latent gene modules. See details.
#' @param M0 Number of clusters. See details.
#' @param clus Initial clustering of scRNASeq data. If NULL, the function will use PCA and Kmeans to do clustering initially.
#' @param K The number of PCs used in the initial clustering. Default = 20.
#' @param iter Number of EM iterations for full data set. See details.
#' @param est_z The iteration starts to update z.
#' @param impt_it The iteration starts to sample new imputed values in initial phase. See details.
#' @param max_lambda Whether to maximize over lambda.
#' @param est_lam The iteration starts to estimate lambda.
#' @param penl L1 penalty for the factor loadings.
#' @param sigma0 The variance of the prior distribution of \eqn{\mu}.
#' @param pi_alpha The hyperparameter of the prior distribution of \eqn{\pi}. See details.
#' @param beta A G by K0 matrix. Initial values for factor loadings (B). If null, beta will initialze from normal distribution with mean zero and variance M0/K0. See details.
#' @param lambda A M0 by K0 matrix. Initial values for the variances of factors. Each column is for a cell cluster. If null, lambda will initialize to be 1/M0. See details.
#' @param sigma A G by M0 matrix. Initial values for the variance of idiosyncratic noises. Each column is for a cell cluster. If null, sigma will initialize to be 1. See details.
#' @param mu A G by M0 matrix. Initial values for the gene expression mean of each cluster. Each column is for a cell cluster. If NULL, it will take the sample mean of cells weighted by the probability in each cluster. See details.
#' @param p_min Initialize parameters using genes expressed in at least \emph{p_min} proportion of cells. If the number genes selected is less than \emph{min_gene}, select \emph{min_gene} genes with higest proportion of non zeros. Default = 0.8.
#' @param min_gene Minimal number of genes used in the initial phase. See details.
#' @param cutoff The value below cutoff is treated as no expression. Default = 0.1.
#' @param verbose Whether to show some intermediate results. Default = False.
#' @param num_mc The number of Gibbs steps to generate new imputed data.
#' @return \code{SIMPLE} returns a list of results in the following order.
#' \enumerate{
#' \item{loglik}{The log-likelihood of the imputed gene expression at each iteration.}
#' \item{pi}{Probabilites of cells belong to each cluster.}
#' \item{mu}{Mean expression for each cluster}
#' \item{sigma}{Variances of idiosyncratic noises for each cluster.}
#' \item{beta}{Factor loadings.}
#' \item{lambda}{Variances of factors for each cluster.}
#' \item{z}{The probability of each cell belonging to each cluster.}
#' \item{Ef}{Conditonal expection the factors for each cluster \eqn{E(f_i|z_i = m)}.
#' A list with length M0, each element in the list is a n by K0 matrix.}
#' \item{Varf}{Conditonal covariance of factors for each cluster \eqn{Var(f_i|z_i = m)}.
#' A list with length M0, each element in the list is a K0 by K0 matrix.}
#' \item{Yimp0}{A matrix contains the expectation of imputed expression.}
#' \item{Y}{Last sample of imputed matrix.}
#' \item{pg}{A G by M0 matrix, dropout rate for each gene in each cluster.}
#' \item{geneM}{Gene mean. If centerized each gene before estimating the parameters, provide the overall mean of gene expression removed from the data matrix. }
#' \item{geneSd}{Gene standard deviation. If scaled each gene before estimating the parameters, provide the overall standard deviation of gene expression removed from the data matrix. }
#' \item{initclus}{Output initial cluster results.}
#' }
#' @seealso [SIMPLE_B()]
#' @examples
#' library(foreach) \cr
#' library(doParallel) \cr
#' library(SIMPLE) \cr
#' source("SIMPLE/utils.R") \cr
#'
#' # simulate number of clusters \cr
#' M0 = 3 \cr
#' # number of cells \cr
#' n = 300 \cr
#' simu_data = simulation_bulk(n=300, S0 = 20, K = 6, MC=M0, block_size = 32, indepG = 1000 - 32*6, verbose=F, overlap=0) \cr
#' Y2 = simu_data$Y2 \cr
#' K0 = 6 # number of factors \cr
#' registerDoParallel(cores = 6)  # parallel \cr
#' # estimate the parameters \cr
#' result <- scimpclu(Y2, K0, M0, celltype=rep(1, n), clus = NULL, K = 20, p_min = 0.5, max_lambda=T, min_gene = 200,cutoff=0.01) \cr
#' # sample imputed values \cr
#' result2 = do_impute(Y2, result$Y, result$beta, result$lambda, result$sigma, result$mu, result$pi, result$geneM, result$geneSd, rep(1, n), mcmc=50, burnin = 5, pg = result$pg, cutoff = 0.01) \cr
#' # evaluate cluster performance \cr
#' celltype_true = simu_data$Z \cr
#' mclust::adjustedRandIndex(apply(result$z,1, which.max), celltype_true) \cr
#' # or redo clustering based on imputed values (sometimes work better for real data) \cr
#' getCluster(result2$impt, celltype_true, Ks = 20, M0 = M0)[[1]] \cr
#'
#' @author Zhirui Hu, \email{zhiruihu@g.harvard.edu}
#' @author Songpeng Zu, \email{songpengzu@g.harvard.edu}
#' @export
SIMPLE <- function(dat, K0, M0 = 1, iter = 10, est_lam = 1, impt_it = 5, penl = 1, sigma0 = 100, pi_alpha = 1, beta = NULL, verbose = F, max_lambda = F, lambda = NULL, sigma = NULL, mu = NULL, est_z = 1, clus = NULL, p_min = 0.8, cutoff = 0.5, K = 10, min_gene = 300, num_mc = 3, fix_num = F, clus_opt = 2, lower = -Inf, upper = Inf) {
  # EM algorithm
  # initiation
  G <- nrow(dat)
  n <- ncol(dat)
  z <- NULL
  Y <- dat # imputed matrix
  # gene_mean = rep(0,G) # remove gene mean for factor, otherwise may be confounded with B, i.e. u = B1, f1=rep(1,n)
  pg <- matrix(p_min, G, M0)

  pi <- rep(1 / M0, M0) # prob of z

  # random start
  if (is.null(lambda)) lambda <- matrix(1 / M0, M0, K0) # sum to M0

  if (is.null(mu)) mu <- matrix(0, G, M0)
  if (is.null(sigma)) sigma <- matrix(1, G, M0)

  if (is.null(beta)) beta <- matrix(rnorm(G * K0), G, K0) / sqrt(K0) * sqrt(M0)


  # # init clustering
  # if(is.null(clus))
  # {
  #   s = svd(t(scale(t(dat))))
  #   km0 <- kmeans(s$v[,1:K], M0, iter.max = 80, nstart = 300)
  #   clus = km0$cluster
  #   if(verbose) {
  #     print(mclust::adjustedRandIndex(clus, celltype_true))
  #     print(xtabs(~clus+celltype_true))
  #   }
  # }

  # inital impution only for low dropout genes
  n1 <- rowMeans(dat > cutoff)
  if (fix_num) {
    hq_ind <- order(n1, decreasing = T)[1:min_gene]
  }
  else {
    hq_ind <- which(n1 >= p_min) #
    if (length(hq_ind) < min_gene) hq_ind <- order(n1, decreasing = T)[1:min_gene] # fix number of hq genes for simulation, need to change back
  }

  print(paste("inital impution for ", length(hq_ind), "high quality genes")) # low dropout

  # init clustering
  if (is.null(clus)) {
    Y2_scale <- t(scale(t(dat[hq_ind, ])))
    s <- svd(Y2_scale)
    if (clus_opt == 1) {
      km0 <- kmeans(t(Y2_scale) %*% s$u[, 1:K], M0, iter.max = 80, nstart = 300) # for high dropout rate
    }
    else {
      km0 <- kmeans(s$v[, 1:K], M0, iter.max = 80, nstart = 300)
    }
    clus <- km0$cluster
    if (verbose & !is.null(celltype_true)) {
      print(mclust::adjustedRandIndex(clus, celltype_true))
      print(xtabs(~ clus + celltype_true))
    }
  }
  z <- matrix(0, n, M0)
  for (m in 1:M0) z[clus == m, m ] <- 1


  if (is.null(clus)) {
    res <- init_impute(dat[hq_ind, ], 1, rep(1, n), p_min, cutoff = cutoff, verbose = verbose) # ???[hq_ind, ]
    res[[2]] <- res[[2]] %*% t(rep(1, M0))
  } else {
    res <- init_impute(dat[hq_ind, ], M0, clus, p_min, cutoff = cutoff, verbose = F) # verbose
  }


  print("impute for hq genes")
  impute_hq <- EM_impute(res[[1]], dat[hq_ind, ], res[[2]], M0, K0, cutoff, 20, beta[hq_ind, ], sigma[hq_ind, , drop = F], lambda, pi, z, mu = NULL, celltype = clus, penl, est_z, max_lambda, est_lam, impt_it, sigma0, pi_alpha, verbose = verbose, num_mc = num_mc, lower = lower, upper = upper) # iter, M0=1?


  pg[hq_ind, ] <- res[[2]]
  beta[hq_ind, ] <- impute_hq$beta
  sigma[hq_ind, ] <- impute_hq$sigma
  mu[hq_ind, ] <- impute_hq$mu # check this, as imputed mu may not be zero
  Y[hq_ind, ] <- impute_hq$Y
  # gene_mean[hq_ind] = impute_hq$geneM
  z <- impute_hq$z


  nz <- colSums(impute_hq$z)
  Vm <- lapply(1:M0, function(m) impute_hq$Varf[[m]] * nz[m])

  # inital beta for other genes
  print("initial estimate beta fpr lq genes:")
  lq_ind <- setdiff(1:G, hq_ind) # which(n1 < p_min )
  # estimate beta and impute: only for positive part? (only impute for genes with more than 10% nonzero)

  # M step also estimate mu
  res <- foreach(g = lq_ind, .combine = rbind) %dopar% { # sapply(lq_ind, function(g){#
    V <- 0
    for (m in 1:M0) V <- V + Vm[[m]] / sigma[g, m]
    W_temp <- c()
    Y_temp <- c()

    for (m in 1:M0)
    {
      Y_temp <- c(Y_temp, dat[g, ] * sqrt(z[, m]) / sqrt(sigma[g, m]))
      Wb <- impute_hq$Ef[[m]] * sqrt(z[, m]) / sqrt(sigma[g, m])
      Wmu <- matrix(0, n, M0)
      Wmu[, m] <- sqrt(z[, m]) / sqrt(sigma[g, m]) # n * M
      W_temp <- rbind(W_temp, cbind(Wmu, Wb))
    }

    ML <- cbind(matrix(0, K0, M0), chol(V))

    W_aug <- rbind(W_temp, ML) # (n+K) * (M + K)
    Y_aug <- c(Y_temp, rep(0, K0)) # G*(n+K)


    penalty <- penl / (M0 * n + K0) / var(Y_aug) # sigma^2
    fit1m <- glmnet(W_aug, Y_aug,
      family = "gaussian", alpha = 1, intercept = F, standardize = F, nlambda = 1, lambda = penalty * K0 / (M0 + K0),
      penalty.factor = c(rep(0, M0), rep(1, K0))
    ) # K dimensional, n+K data

    coeff <- fit1m$beta[-1:-M0, 1]
    tempmu <- fit1m$beta[1:M0, 1]

    # sg =  sum((Y_temp - fit1m$a0 - coeff %*% t(W_temp))^2) + sum(( coeff %*%V)*  coeff)* ns
    #
    # c(fit1m$a0,  coeff, (sg+1)/(ns + 3))

    sg <- sapply(1:M0, function(m) {
      (sum((dat[g, ] - tempmu[m] - coeff %*% t(impute_hq$Ef[[m]]))^2 * z[, m]) + sum((coeff %*% Vm[[m]]) * coeff))
    })
    c(fit1m$beta[, 1], rep((sum(sg) + 1) / (n + 3), M0))
  }

  mu[lq_ind, ] <- matrix(res[, 1:M0], ncol = M0)
  beta[lq_ind, ] <- matrix(res[, (M0 + 1):(K0 + M0)], ncol = K0)

  sigma[lq_ind, ] <- matrix(res[, -1:-(K0 + M0)], ncol = M0)
  sigma[sigma > 9] <- 9
  sigma[sigma < 1e-4] <- 1e-4


  # imputation set dropout rate as p_min
  if (M0 > 1) {
    im <- apply(z, 1, function(x) which(rmultinom(1, 1, x) == 1)) # sample membership
  } else {
    im <- rep(1, n)
  }
  for (i in 1:n)
  {
    m <- im[i]
    # vr = impute_hq$Varf[[m]]
    # f_i = rmvnorm(1, impute_hq$Ef[[m]][i,], vr)

    ind <- which(dat[lq_ind, i] <= cutoff)
    ind <- lq_ind[ind]

    ms <- mu[ind, m] + beta[ind, , drop = F] %*% impute_hq$Ef[[m]][i, ]
    sds <- sqrt(sigma[ind, m])
    p <- pg[ind, clus[i]] # need celltype

    prob <- pnorm(cutoff, mean = ms, sd = sds) # compute x<0 prob
    prob_drop <- (1 - p) / (prob * p + (1 - p))
    I_drop <- rbinom(length(ind), 1, prob_drop)


    # imputation for dropout
    impt <- rep(0, length(ind))
    impt[I_drop == 1] <- rnorm(sum(I_drop == 1), ms[I_drop == 1], sds[I_drop == 1])

    # imputation for non-dropout
    if (sum(I_drop == 0) > 0) {
      impt[I_drop == 0] <- rtnorm(sum(I_drop == 0), upper = cutoff, mean = ms[I_drop == 0], sd = sds[I_drop == 0])
    }

    Y[ind, i] <- impt
  }


  print("impute for all genes")
  impute_result <- EM_impute(Y, dat, pg, M0, K0, cutoff, iter, beta, sigma, impute_hq$lambda, impute_hq$pi, impute_hq$z, mu = NULL, celltype = clus, penl, est_z, max_lambda, est_lam, impt_it = 1, sigma0, pi_alpha, verbose = verbose, num_mc = num_mc, lower = lower, upper = upper)


  impute <- matrix(0, n, G)
  for (m in 1:M0)
  {
    impute <- impute + t(impute_result$mu[, m] + impute_result$beta %*% t(impute_result$Ef[[m]])) * impute_result$z[, m]
  }
  impute <- t(impute) * impute_result$geneSd + impute_result$geneM

  return(list("loglik" = impute_result$loglik, "pi" = impute_result$pi, "mu" = impute_result$mu, "sigma" = impute_result$sigma, "beta" = impute_result$beta, "lambda" = impute_result$lambda, "z" = impute_result$z, "Ef" = impute_result$Ef, "Varf" = impute_result$Varf, "geneM" = impute_result$geneM, "geneSd" = impute_result$geneSd, "Yimp0" = impute, "Y" = impute_result$Y, "pg" = pg, "initclus" = clus))
}
