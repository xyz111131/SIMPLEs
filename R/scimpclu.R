#Factor mixture model, laplace prior for B
# full data
# library(ggplot2)
# library(reshape2)
# library(gridExtra)
# library(MASS)
# library(glmnet)
# library(matrixStats)
# library(mixtools)
# library(Rsolnp)
# library(msm)

PI = pi

#source('R/fit_ztrunc_0228.R')
#source('R/EM_impute.R')

# imputation
init_impute <- function(Y2,M0, clus,p_min =0.6, cutoff = 0.5,verbose=F)
{
  # fit truncnorm for each cluster
  impute=Y2
  G = nrow(Y2)
  pg = matrix(0,G, M0)
  for(i in 1:M0)
  {
    temp_dat = as.matrix(Y2[,clus==i])
    result <- ztruncnorm(temp_dat, cutoff = cutoff, p_min = p_min) #, s_upper=3
    mu = result[[1]][,1]
    sd = result[[1]][,2]
    pg[,i] = result[[2]]

    pg[pg > 0.99] = 0.99

    ind = which(is.na(mu))
    if(verbose) print(paste(length(ind), "not fit TN"))
    mu[ind] = rowMeans(Y2[ind,clus==i, drop=F])
    sd[ind] = rowSds(Y2[ind,clus==i, drop=F])
    pg[ind, i] = 1

    if(verbose) print(summary(sd))

    for(j in which(clus==i))
    {
      ind = which(Y2[,j] <=cutoff)
      # impute Y2
      ms = mu[ind]
      sds = sd[ind]
      p = pg[ind, i]
      prob = pnorm(cutoff, mean = ms, sd = sds) # compute x<0 prob
      prob_drop = (1-p) / (prob * p + (1-p))
      I_drop = rbinom(length(prob_drop), 1, prob_drop)

      # imputation for dropout
      impute[ind[I_drop==1],j] = rnorm(sum(I_drop==1), ms[I_drop==1], sds[I_drop==1])
      #åimpute[ind[I_drop==1],j] = ms[I_drop==1]

      # imputation for non-dropout
      if(sum(I_drop==0) >0) {
        #r = (cutoff - ms[I_drop==0])/sds[I_drop==0]
        impute[ind[I_drop==0],j] = rtnorm(sum(I_drop==0), upper = cutoff, mean = ms[I_drop==0], sd = sds[I_drop==0]) #ms[I_drop==0] - sds[I_drop==0] * dnorm(r)/pnorm(r)
      }
    }
  }

  impute[is.na(impute)] = 0

  #print(sum(impute < -2))
  if(verbose) print(summary(c(impute[Y2==0])))
  #impute[impute < -3] = -3

  return(list(impute, pg))
}

#' @export
scimpclu <- function(dat, K0, M0 = 1, iter=10, est_lam = 2, impt_it = 5, penl =1 ,sigma0 = 100,  pi_alpha = 1, beta = NULL, verbose=F, max_lambda=F, lambda =NULL, sigma=NULL, mu=NULL, est_z = 1, clus = NULL,p_min = 0.8, cutoff=0.5, K = 10, min_gene = 300, num_mc = 3, lower = -Inf, upper = Inf)  # M0 is only for inital impute for hq genes, celltype = NULL,
{
  #EM algorithm
  #initiation
  G = nrow(dat)
  n= ncol(dat)
  z = NULL
  Y = dat # imputed matrix
  #gene_mean = rep(0,G) # remove gene mean for factor, otherwise may be confounded with B, i.e. u = B1, f1=rep(1,n)
  pg = matrix(p_min, G, M0)

  pi = rep(1/M0, M0) # prob of z

  # random start
  if(is.null(lambda)) lambda =matrix(1/M0, M0, K0) #sum to M0

  if(is.null(mu)) mu =  matrix(0, G, M0)
  if(is.null(sigma)) sigma <-  matrix(1, G, M0)

  if(is.null(beta)) beta <- matrix(rnorm(G*K0), G, K0)/sqrt(K0)*sqrt(M0)


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
  n1 = rowMeans(dat > cutoff)
  hq_ind = which(n1 >= p_min ) #
  if(length(hq_ind) < min_gene) hq_ind = order(n1, decreasing = T)[1:min_gene]
  print(paste("inital impution for ", length(hq_ind), "high quality genes" )) #low dropout

  # init clustering
  if(is.null(clus))
  {
    Y2_scale = t(scale(t(dat[hq_ind, ])))
    s = svd(Y2_scale)
    km0 <- kmeans(t(Y2_scale)%*% s$u[,1:K], M0, iter.max = 80, nstart = 300) #s$v[,1:K]
    #km0 <- kmeans(s$v[,1:K], M0, iter.max = 80, nstart = 300)
    clus = km0$cluster
    if(verbose & !is.null(celltype_true)) {
      print(mclust::adjustedRandIndex(clus, celltype_true))
      print(xtabs(~clus+celltype_true))
    }
  }
  z = matrix(0, n, M0)
  for(m in 1:M0) z[clus==m, m ] =1


  if(is.null(clus))
  {
    res = init_impute(dat[hq_ind, ],1, rep(1,n), p_min, cutoff = cutoff, verbose = verbose) # ???[hq_ind, ]
    res[[2]] = res[[2]] %*% t(rep(1, M0))
   }else{
     res = init_impute(dat[hq_ind, ] , M0, clus, p_min, cutoff = cutoff, verbose = F) #verbose
   }


  print("impute for hq genes")
  impute_hq <- EM_impute(res[[1]], dat[hq_ind, ], res[[2]], M0, K0, cutoff, 20, beta[hq_ind,], sigma[hq_ind,,drop=F], lambda, pi, z, mu = NULL, celltype = clus, penl, est_z, max_lambda, est_lam, impt_it, sigma0, pi_alpha,verbose=verbose, num_mc = num_mc, lower = lower, upper = upper) #iter, M0=1?


  pg[hq_ind, ] = res[[2]]
  beta[hq_ind, ] = impute_hq$beta
  sigma[hq_ind, ]= impute_hq$sigma
  mu[hq_ind, ] = impute_hq$mu  # check this, as imputed mu may not be zero
  Y[hq_ind, ] = impute_hq$Y
  #gene_mean[hq_ind] = impute_hq$geneM
  z = impute_hq$z


  nz =colSums(impute_hq$z)
  Vm = lapply(1:M0, function(m)  impute_hq$Varf[[m]] * nz[m])

  # inital beta for other genes
  print("initial estimate beta fpr lq genes:")
  lq_ind = setdiff(1:G, hq_ind) #which(n1 < p_min )
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


  # imputation set dropout rate as p_min
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
    p = pg[ind, clus[i]]  # need celltype

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



  # if(verbose){
  #   par(mfrow=c(2,2), pch= 16)
  #   plot(dat[hq_ind[5],], Y[hq_ind[5],], xlab = "data", ylab = "imputed", main = "hq gene"); abline(c(0,1), col=2)
  #   plot(dat[lq_ind[5],], Y[lq_ind[5],], xlab = "data", ylab = "imputed", main = "lq gene"); abline(c(0,1), col=2)
  #
  #   k = which.max(abs(beta[lq_ind[5],]))
  #   ef = sapply(1:n, function(i) impute_hq$Ef[[im[i]]][i,k])
  #   plot(ef, Y[lq_ind[5],], col=celltype_true, ylab="imputed", xlab="F"); abline(h=mu[lq_ind[5],], col=1:M0, lty=2)
  #   plot(ef, dat[lq_ind[5],],col=celltype_true, ylab="data", xlab="F")
  #
  #
  #   hist(sigma[,1], xlab = "sigma")
  #   # boxplot(Y[5, ]~celltype_true)
  #   # boxplot(dat[5, ]~celltype_true)
  # }


  print("impute for all genes")
  impute_result <- EM_impute(Y, dat, pg, M0, K0, cutoff, iter, beta, sigma, impute_hq$lambda, impute_hq$pi, impute_hq$z, mu=NULL, celltype = clus, penl, est_z, max_lambda, est_lam, impt_it=1, sigma0, pi_alpha,verbose=verbose,num_mc = num_mc, lower = lower, upper = upper)


  impute = matrix(0, n, G)
  for(m in 1:M0)
  {
    impute = impute + t(impute_result$mu[,m] + impute_result$beta %*% t(impute_result$Ef[[m]])) *impute_result$z[,m]
  }
  impute = t(impute)*impute_result$geneSd + impute_result$geneM

  return(list("loglik" = impute_result$loglik, "pi" = impute_result$pi, "mu" = impute_result$mu, "sigma" = impute_result$sigma, "beta" = impute_result$beta, "lambda"= impute_result$lambda, "z" = impute_result$z, "Ef" = impute_result$Ef, 'Varf' = impute_result$Varf, 'geneM' = impute_result$geneM, 'geneSd'=impute_result$geneSd, 'Yimp0' =impute, 'Y' = impute_result$Y, 'pg' = pg, 'initclus' = clus ))
}


