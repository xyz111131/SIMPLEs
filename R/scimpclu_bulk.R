# Factor mixture model, laplace prior for B
# integrate bulk

#' @import MASS
#' @import glmnet
#' @import matrixStats
#' @import Rsolnp
#' @importFrom mixtools rmvnorm
#' @importFrom msm rtnorm

PI = pi

#source('R/fit_ztrunc_0228.R')
#source('R/EM_impute.R')

init_impute_bulk <- function(Y2, clus, bulk, pg1, cutoff = 0.1,verbose=F)
{
  # fit truncnorm for each cluster
  impute0 =Y2
  G = nrow(Y2)
  pg = pg1 #matrix(0,G, M0)
  #ndropout = pg1
  M0 = max(clus)
  for(i in 1:M0)
  {
    temp_dat = as.matrix(Y2[,clus==i])
    #p0 = rowMeans(temp_dat <= cutoff)
    impute = temp_dat

    result <- ztruncnorm(temp_dat, cutoff = cutoff, p_min = 0.01, p_max = pg1[,i]) #pg1 >= 0.1
    mu = result[[1]][,1]
    sd = result[[1]][,2]
    pg[,i] = result[[2]]


    ind = which(is.na(mu))
    if(verbose) print(paste(length(ind), "not fit TN"))
    mu[ind] = rowMeans(Y2[ind,clus==i, drop=F])
    sd[ind] = rowSds(Y2[ind,clus==i, drop=F])
    pg[ind, i] =pg1[ind, i] #1


    if(verbose) {
      par(mfrow=c(2,2), pch=16)
      plot(bulk[,i], mu, col=rgb(1,0,0,0.5), xlab = "bulk", ylab = "est_mu", main= paste("cluster", i))
      a = coef(lm(mu~bulk[,i], weights = bulk[,i]^2))
      abline(a, col=2)

      plot(bulk[,i], rowMeans(temp_dat), col=rgb(0,1,0,0.5), xlab = "bulk", ylab = "raw_mean")
      b = coef(lm(rowMeans(temp_dat)~bulk[,i], weights = bulk[,i]^2))
      abline(b, col=2)

      # dropout rate
      plot(pg1[,i], pg[,i], col=rgb(0,1,1,0.5), xlab = "amplified rate", ylab = "prob normal component")
    }


    # imputation for each gene
    for(j in 1:G)
    {
      imp = which(temp_dat[j,] <=cutoff)
      if(length(imp) == 0) next;
      prob = pnorm(cutoff, mean = mu[j], sd = sd[j]) # compute x<0 prob

      # imputation for dropout from norm, zero comp, censored part
      prob_drop = (pg[j,i] / pg1[j,i] -  pg[j,i]) / (prob * pg[j,i] + (1-pg[j,i])) # q * d, pg1 = 1-d, p = (1-d)*q, q = pg/pg1
      prob_zero = (1 - pg[j,i] / pg1[j,i]) / (prob * pg[j,i] + (1-pg[j,i])) #1-q
      prob_trunc = prob * pg[j,i] / (prob * pg[j,i] + (1-pg[j,i])) #q(1-d)*prob

      I_drop = rmultinom(length(imp), 1, c(prob_drop, prob_zero, prob_trunc))

      # imputation for dropout
      impute[j,imp[I_drop[1,] == 1]] = rnorm(sum(I_drop[1,]==1), mu[j], sd[j])

      # imputation for others
      impute[j,imp[I_drop[2,] == 1]] = rnorm(sum(I_drop[2,]==1), -1, 0.5)

      impute[j, imp[I_drop[3,]==1]] = rtnorm(sum(I_drop[3,]==1), upper = cutoff, mean = mu[j], sd=sd[j])
    }

    #print(impute[j,imp])


    impute0[, clus==i] = impute;

    if(verbose) {
      plot(bulk[,i], rowMeans(impute), col=rgb(0,0,1,0.5), xlab = "bulk", ylab = "imputed mean");
      abline(b, col=2);

    }
  }

  impute0[is.na(impute0)] = 0

  ##plot(Y2[803 ,clus==1], impute0[803 ,clus==1])

  return(impute0) #ndropout, pg
}

#' @export
scimpclu_bulk <- function(dat, K0, bulk, M0=1, celltype = NULL, clus = NULL, K = 20, iter=10, est_z = 1,  impt_it = 5, max_lambda=F, est_lam = 1,penl =1 ,sigma0 = 100,  pi_alpha = 1, beta = NULL, lambda =NULL, sigma=NULL, mu=NULL, p_min = 0.8, min_gene = 300, cutoff=0.1, verbose=F, num_mc = 3, fix_num = F, clus_opt = 2)  # iter: EM step for all genes
{
  #EM algorithm
  #initiation
  G = nrow(dat)
  n= ncol(dat)
  Y = dat # imputed matrix
  gene_mean = rep(0,G) # remove gene mean for factor, otherwise may be confounded with B, i.e. u = B1, f1=rep(1,n)

  pi = rep(1/M0, M0) # prob of z

  # random start
  if(is.null(lambda)) lambda = matrix(1/M0, M0, K0) #sum M0

  if(is.null(mu)) mu =  matrix(0, G, M0)
  if(is.null(sigma)) sigma <-  matrix(1, G, M0)

  if(is.null(beta)) beta <- matrix(rnorm(G*K0), G, K0)/sqrt(K0) *sqrt(M0)


  # inital impution only for low dropout genes
  n1 = rowMeans(dat > cutoff)
  if(fix_num)
  {
    hq_ind = order(n1, decreasing = T)[1:min_gene]
  }else{
    hq_ind = which(n1 >= p_min)
    if(length(hq_ind) < min_gene) hq_ind = order(n1, decreasing = T)[1:min_gene] # need to change back
  }
  # regression for hq genes and estimate dropout rate for all genes
  MB = max(celltype)
  pg = matrix(0, G, MB)
  for(i in 1:MB)
  {
    #b = coef(lm(rowMeans(dat[hq_ind, celltype==i])~bulk[hq_ind,i], weights = bulk[hq_ind,i]^2))
    b = c(0,1)
    pg[,i] = (rowMeans(dat[,celltype==i ]) + p_min) / (1 + b[1] + b[2] * bulk[,i])
  }

  pg[pg>1] = 1
  pg[pg<0.1] = 0.1

  print(paste("initial impution for hq genes: ", length(hq_ind))) #low dropout

  res = init_impute_bulk(dat[hq_ind, ], celltype, bulk[hq_ind,,drop=F ], pg[hq_ind, ,drop=F], cutoff = cutoff, verbose = F) #verbose

  # init clustering
  if(is.null(clus))
  {
    #s = svd(t(scale(t(res))))
    #km0 <- kmeans(s$v[,1:K], M0, iter.max = 80, nstart = 300)

    Y2_scale = t(scale(t(res)))
    s = svd(Y2_scale)
    if(clus_opt==1)
    {
      km0 <- kmeans(t(Y2_scale)%*% s$u[,1:K], M0, iter.max = 80, nstart = 300)
    }else{
      km0 <- kmeans(s$v[,1:K], M0, iter.max = 80, nstart = 300)
    }
    clus = km0$cluster

    if(verbose) {
      print(xtabs(~clus))
      print("initial clustering randInd: ")
      if(!is.null(celltype))
      {
        print(mclust::adjustedRandIndex(clus, celltype))
      }
      if(!is.null(celltype_true))
      {
        print(mclust::adjustedRandIndex(clus, celltype_true))
      }

    }

  }

  z = matrix(0, n, M0)
  for(m in 1:M0) z[clus==m, m ] =1

  print("initial estimate factors: ")
  impute_hq <- EM_impute(res, dat[hq_ind, ], pg[hq_ind, , drop=F] , M0, K0, cutoff, 30, beta[hq_ind,], sigma[hq_ind,,drop=F], lambda, pi, z, mu = NULL, celltype = celltype, penl, est_z, max_lambda, est_lam, impt_it, sigma0, pi_alpha,verbose=verbose, num_mc = num_mc)  # iter = 30, mu = NULL


  beta[hq_ind, ] = impute_hq$beta
  sigma[hq_ind, ]= impute_hq$sigma
  mu[hq_ind, ] = impute_hq$mu  # check this, as imputed mu may not be zero
  Y[hq_ind, ] = impute_hq$Y
  gene_mean[hq_ind] = impute_hq$geneM
  z = impute_hq$z


  nz =colSums(impute_hq$z)
  Vm = lapply(1:M0, function(m)  impute_hq$Varf[[m]] * nz[m])

  # inital beta for other genes
  print("initial estimate beta fpr lq genes:")
  #V = impute_hq$Varf[[1]]
  #ML <- chol(V)
  #lq_ind = which(n1 < p_min )
  lq_ind = setdiff(1:G, hq_ind)
  # estimate beta and impute: only for positive part? (only impute for genes with more than 10% nonzero)

  # M step also estimate mu
  res = foreach(g=lq_ind, .combine = rbind) %dopar% {#sapply(lq_ind, function(g){#
    V = 0
    for(m in 1:M0) V = V + Vm[[m]]/sigma[g,m]
    W_temp = c()
    Y_temp = c()

    for(m in 1:M0)
    {
      Y_temp = c(Y_temp, dat[g, ]*sqrt(z[,m])/sqrt(sigma[g,m]))
      Wb = impute_hq$Ef[[m]] *sqrt(z[,m])/sqrt(sigma[g,m])
      Wmu = matrix(0, n, M0)
      Wmu[, m] = sqrt(z[,m])/sqrt(sigma[g,m]) # n * M
      W_temp = rbind(W_temp, cbind(Wmu, Wb))
    }

    ML <- cbind(matrix(0, K0, M0), chol(V))

    W_aug <- rbind(W_temp, ML) #(n+K) * (M + K)
    Y_aug <- c(Y_temp, rep(0,K0)) #G*(n+K)



    # #s_ind = which(dat[g, ] > cutoff)
    # W_temp = impute_hq$Ef[[1]][, ,drop=F] #s_ind
    # Y_temp =  dat[g, ] #res[[1]][g, ], s_ind

    # ns = n #length(s_ind)
    # W_aug <- rbind(W_temp, ML*sqrt(ns)) #(n+K) * K
    # Y_aug <- c(Y_temp, rep(0,K0)) #G*(n+K)

    penalty <- penl/(M0*n + K0)/var(Y_aug) # sigma^2
    fit1m=glmnet(W_aug,Y_aug,family="gaussian",alpha=1, intercept=F, standardize = F, nlambda =1, lambda = penalty * K0/(M0 + K0),
                 penalty.factor = c(rep(0, M0), rep(1, K0))) #K dimensional, n+K data

    coeff = fit1m$beta[-1:-M0,1]
    tempmu = fit1m$beta[1:M0,1]

    # sg =  sum((Y_temp - fit1m$a0 - coeff %*% t(W_temp))^2) + sum(( coeff %*%V)*  coeff)* ns
    #
    # c(fit1m$a0,  coeff, (sg+1)/(ns + 3))

    sg = sapply(1:M0, function(m) {
      (sum( (dat[g, ] - tempmu[m] - coeff%*% t(impute_hq$Ef[[m]]) )^2 *z[,m] ) + sum((coeff %*%Vm[[m]])* coeff))
    })
    c(fit1m$beta[,1], rep((sum(sg)+1)/(n + 3),M0))
  }

  mu[lq_ind, ] = matrix(res[, 1:M0], ncol=M0)
  beta[lq_ind, ] = matrix(res[,(M0+1):(K0+M0)], ncol=K0)

  sigma[lq_ind, ] = matrix(res[, -1:-(K0+M0)], ncol=M0)
  sigma[sigma>9]=9
  sigma[sigma<1e-4]=1e-4


  # imputation set dropout rate as pg
  if(M0 > 1) {
    im =  apply(z,1, function(x) which(rmultinom(1, 1, x)==1)) # sample membership
  }else{
    im = rep(1,n)
  }
  for(i in 1:n)
  {
    m = im[i]
    # vr = impute_hq$Varf[[m]]
    # f_i = rmvnorm(1, impute_hq$Ef[[m]][i,], vr)

    ind = which(dat[lq_ind,i] <=cutoff)
    ind = lq_ind[ind]

    ms = mu[ind,m] + beta[ind,, drop=F] %*% impute_hq$Ef[[m]][i,]
    sds = sqrt(sigma[ind,m])
    p = pg[ind, celltype[i]]  # need celltype

    prob = pnorm(cutoff, mean = ms, sd = sds) # compute x<0 prob
    prob_drop = (1-p) / (prob * p + (1-p))
    I_drop = rbinom(length(ind), 1, prob_drop)


    # imputation for dropout
    impt = rep(0, length(ind))
    impt[I_drop==1] = rnorm(sum(I_drop==1), ms[I_drop==1], sds[I_drop==1])

    # imputation for non-dropout
    if(sum(I_drop==0) >0) {
      impt[I_drop==0] = rtnorm(sum(I_drop==0), upper = cutoff, mean = ms[I_drop==0], sd = sds[I_drop==0])
    }

    Y[ind, i] = impt
  }


  impute = matrix(0, n, G)
  if(iter > 0)
  {
    # EM for all genes
    #Y[hq_ind, ] = Y[hq_ind, ] + impute_hq$geneM
    print("impute for all genes")
    impute_result <- EM_impute(Y, dat, pg, M0, K0, cutoff, iter, beta, sigma, impute_hq$lambda, impute_hq$pi, impute_hq$z, mu=NULL, celltype = celltype, penl, est_z, max_lambda, est_lam, impt_it=1, sigma0, pi_alpha, verbose=verbose, num_mc = num_mc)


    for(m in 1:M0)
    {
      impute = impute + t(impute_result$mu[,m] + impute_result$beta %*% t(impute_result$Ef[[m]])) *impute_result$z[,m]
    }

    impute = t(impute)*impute_result$geneSd + impute_result$geneM

    return(list("loglik" = impute_result$loglik, "pi" = impute_result$pi, "mu" = impute_result$mu, "sigma" = impute_result$sigma, "beta" = impute_result$beta, "lambda"= impute_result$lambda, "z" = impute_result$z, "Ef" = impute_result$Ef, 'Varf' = impute_result$Varf, 'Yimp0' =impute, 'Y' = impute_result$Y, 'pg' = pg,'geneM' = impute_result$geneM, 'geneSd'=impute_result$geneSd))

  }else{
    #gene_mean[lq_ind] = mu[lq_ind,]
    #Y[lq_ind, ] = Y[lq_ind, ] - mu[lq_ind]
    #mu[lq_ind] = 0


    for(m in 1:M0)
    {
      impute = impute + t(mu[,m] + beta %*% t(impute_hq$Ef[[m]])) *impute_hq$z[,m]
    }

    impute = t(impute) + gene_mean

    return(list("loglik" = impute_hq$loglik, "pi" = pi, "mu" = mu, "sigma" = sigma, "beta" = beta, "lambda"= lambda, "z" = impute_hq$z, "Ef" = impute_hq$Ef, 'Varf' = impute_hq$Varf, 'Yimp0' =impute, 'Y' = Y, 'pg' = pg, 'geneM' = gene_mean, 'geneSd'=rep(1,G)))

  }


}




