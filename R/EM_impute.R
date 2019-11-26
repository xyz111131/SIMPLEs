EM_impute <- function(Y, Y0, pg, M0, K0, cutoff, iter, beta, sigma, lambda, pi, z, mu = NULL, celltype = NULL, penl = 1, est_z = 2, max_lambda = T, est_lam = 2, impt_it = 5, sigma0 = 100, pi_alpha = 1, verbose = F, num_mc = 3, lower = -Inf, upper = Inf) # if est_z = 1, not use initial values of z
{
  Y[Y > upper] <- upper
  Y[Y < lower] <- lower

  n <- ncol(Y)
  G <- nrow(Y)

  gene_mean <- rowMeans(Y)
  gene_sd <- rep(1, G) # rowSds(Y)
  Y <- Y - rowMeans(Y) # t(scale(t(Y)))


  if (is.null(mu)) {
    if (M0 > 1 & !is.null(z)) {
      # clus = apply(z, 1, which.max)
      # for(m in 1:M0)
      # {
      #   mu[,m] = rowMeans(Y[,clus==m])
      # }
      mu <- Y %*% z / n
    }
    else {
      mu <- matrix(0, G, M0)
    }
  }

  if (is.null(celltype)) celltype <- rep(1, n) # cell type label for bulk and pg

  loglik <- matrix(0, n, M0)
  loglik0 <- matrix(0, n, M0)
  ## rownames(loglik0) = colnames(dat)
  tot <- rep(0, iter)

  M <- list() # (BSigma^-1B + Lambda^-1)^-1
  Wm <- list() # mean of each factor, K0*n

  equal <- function(x, Em, nz) {
    sum(x)
  }

  fn <- function(x, Em, nz) {
    sum(Em / x + log(x) * (nz - 2)) # -2(alpha-1)
  }

  for (it in 1:iter)
  {
    # E step
    # get expectation of Z
    for (m in 1:M0)
    {
      beta2 <- beta / sigma[, m]
      M[[m]] <- solve(t(beta) %*% beta2 + diag(1 / lambda[m, ]))
      dt <- sum(log(sigma[, m])) + sum(log(lambda[m, ])) - log(det(M[[m]])) # !!! logdet !=det(log=T)
      loglik0[, m] <- apply(Y, 2, function(x) {
        x2 <- x - mu[, m]
        y <- sum(x2 / sigma[, m] * x2)
        x2 <- x2 %*% beta2
        y <- y - sum((x2 %*% M[[m]]) * x2)

        return(y)
      })

      loglik[, m] <- (-loglik0[, m] - dt - G * log(2 * PI)) / 2
    }

    loglik <- t(t(loglik) + log(pi))
    sconst <- rowMaxs(loglik)
    loglik2 <- loglik - sconst
    lik <- exp(loglik2)

    # compute total loglik
    tot[it] <- sum(log(rowSums(lik)) + sconst)

    if (it >= est_z | is.null(z)) {
      z <- lik / rowSums(lik) # n * M0
    }

    nz <- colSums(z)

    # get expectation of W
    for (m in 1:M0)
    {
      Wm[[m]] <- t(M[[m]] %*% t(beta) %*% ((Y - mu[, m]) / sigma[, m])) # n * K0
    }


    # imputation
    if (it >= impt_it + 1) {
      for (iii in 1:num_mc)
      {
        if (M0 > 1) {
          im <- apply(z, 1, function(x) which(rmultinom(1, 1, x) == 1)) # sample membership
        } else {
          im <- rep(1, n)
        }
        for (i in 1:n)
        {
          m <- im[i]
          vr <- M[[m]]
          f_i <- rmvnorm(1, Wm[[m]][i, ], vr)

          ind <- which(Y0[, i] <= cutoff)

          ms <- (mu[ind, m] + beta[ind, ] %*% t(f_i)) * gene_sd[ind] + gene_mean[ind]
          sds <- sqrt(sigma[ind, m]) * gene_sd[ind]
          p <- pg[ind, celltype[i]] # need celltype

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

          impt[impt > upper] <- upper
          impt[impt < lower] <- lower

          Y[ind, i] <- (impt - gene_mean[ind]) / gene_sd[ind]
        }

        # get expectation of Z
        for (m in 1:M0)
        {
          loglik0[, m] <- apply(Y, 2, function(x) {
            x2 <- x - mu[, m]
            y <- sum(x2 / sigma[, m] * x2)
            x2 <- x2 %*% beta2
            y <- y - sum((x2 %*% M[[m]]) * x2)

            return(y)
          })

          loglik[, m] <- (-loglik0[, m] - dt - G * log(2 * PI)) / 2
        }

        loglik <- t(t(loglik) + log(pi))
        sconst <- rowMaxs(loglik)
        loglik2 <- loglik - sconst
        lik <- exp(loglik2)

        if (it >= est_z | is.null(z)) {
          z <- lik / rowSums(lik) # n * M0
        }

        nz <- colSums(z)

        # get expectation of W
        for (m in 1:M0)
        {
          Wm[[m]] <- t(M[[m]] %*% t(beta) %*% ((Y - mu[, m]) / sigma[, m])) # n * K0
        }
      }
    }


    Vm <- lapply(1:M0, function(m) M[[m]] * nz[m])


    # M step
    res <- foreach(g = 1:G, .combine = rbind) %dopar% {
      V <- 0
      for (m in 1:M0) V <- V + Vm[[m]] / sigma[g, m]
      W_temp <- c()
      Y_temp <- c()

      for (m in 1:M0)
      {
        Y_temp <- c(Y_temp, (Y[g, ] - mu[g, m]) * sqrt(z[, m]) / sqrt(sigma[g, m]))
        W_temp <- rbind(W_temp, Wm[[m]] * sqrt(z[, m]) / sqrt(sigma[g, m]))
      }

      ML <- chol(V)

      W_aug <- rbind(W_temp, ML) # (n+K) * K
      Y_aug <- c(Y_temp, rep(0, K0)) # G*(n+K)


      penalty <- penl # sigma^2
      fit1m <- glmnet(W_aug, Y_aug,
        family = "gaussian", alpha = 1, intercept = F, standardize = F, nlambda = 1, lambda = mean(penalty) / (M0 * n + K0) / var(Y_aug),
        penalty.factor = rep(1, K0)
      ) # K dimensional, n+K data

      nb <- fit1m$beta[, 1]
      # max sigma
      # sg = sapply(1:M0, function(m) {
      #  (sum( (Y[g, ] - mu[g,m] - nb%*% t(Wm[[m]]) )^2 *z[,m] ) + sum((nb %*%Vm[[m]])* nb) + 1) / ( nz[m] + 3)
      # }))

      sg <- sapply(1:M0, function(m) {
        (sum((Y[g, ] - mu[g, m] - nb %*% t(Wm[[m]]))^2 * z[, m]) + sum((nb %*% Vm[[m]]) * nb))
      })
      c(nb, rep((sum(sg) + 1) / (n + 3), M0))
    }

    beta <- matrix(res[, 1:K0], ncol = K0)
    sigma <- matrix(res[, -1:-K0], ncol = M0)
    sigma[sigma > 9] <- 9
    sigma[sigma < 1e-4] <- 1e-4

    # max lambda
    if (max_lambda & it >= est_lam) {
      ind_m <- which((nz > K0 + 1) & (nz > 3)) # at least K0 samples to estimate lambda
      mt <- length(ind_m)
      if (mt > 1) {
        Em <- sapply(ind_m, function(m) { # 1:M0
          # a = Wm[[m]]^2 * z0[,m] # windsor, n * K0
          # b = colQuantiles(a[z0[,m] > 0.5, ], probs = 0.95)
          # for(k in 1:K0)
          # {
          #   a[ a[,k] > b[k],k] = b[k]
          #   #a[ a[,k] > b[k,2],k] = b[k,2]
          # }
          # colSums(a) + diag(M[[m]]) * nz0[m]
          colSums(Wm[[m]]^2 * z[, m]) + diag(M[[m]]) * nz[m]
        })
        # print(round(Em,4))

        lambda[-ind_m, ] <- 1 / M0
        lam_init <- t(Em + 1) / (nz[ind_m] + 1)
        lam_init <- t(t(lam_init) / colSums(lam_init)) * mt / M0

        for (k in 1:K0)
        {
          if (all(beta[, k] == 0)) {
            lambda[, k] <- 1 / M0
          } else {
            lambda[ind_m, k] <- solnp(lam_init[, k], fun = fn, eqfun = equal, eqB = mt / M0, LB = rep(0, mt), UB = rep(mt / M0, mt), control = list(inner.iter = 500, trace = 0), nz = nz[ind_m], Em = Em[k, ])$pars
          }
        }
      } else {
        lambda <- matrix(1 / M0, M0, K0)
      }
    }

    # maximizme mu
    for (m in 1:M0)
    {
      sigma_temp <- (nz[m] + sigma[, m] / sigma0)
      beta2 <- beta / sigma_temp

      bigM <- diag(1 / sigma_temp) - beta2 %*% solve(diag(sigma0 / lambda[m, ]) + t(beta2) %*% beta) %*% t(beta2)

      mu[, m] <- bigM %*% (Y %*% z[, m])
    }

    # max pi
    pi <- (nz + pi_alpha - 1) / (n + pi_alpha * M0 - M0)


    # plot((dat[5,]-pos_mean[5])/pos_sd[5], Y[5,]);abline(c(0,1),col=2)
    # boxplot(Y[803, ] + gene_mean[803]~celltype)
    # print(beta[803,])
    # print(sigma[803,])

    # if(impt_it == 1)
    # {
    #   boxplot(Y[1344, ] + gene_mean[1344]~celltype)
    #   print(beta[1344,])
    #   print(sigma[1344,])
    # }




    # plot beta every 10 iter
    if (verbose & it %% 10 == 0) { # == iter
      print(paste(it, tot[it]))
      print("number of cells for each cluster: ")
      print(round(nz))
      print("sum of factor variances for each cluster: ")
      print(round(rowSums(lambda), 2))

      if (!is.null(celltype_true)) {
        # matplot(apply(z, 1, which.max), pch=16, cex=0.5)# cell is ordered by true cell type
        # y = cumsum(xtabs(~celltype_true))
        # for(yt in y) abline(v = yt, col=2, lty=2)
        print("clustering randInd: ")
        print(mclust::adjustedRandIndex(im, celltype_true))
      }

      ## plot(ggplot(melt(beta), aes(Var2,Var1)) +geom_tile(aes(fill = value)) +
      ##       scale_fill_gradient2()+theme_bw()+xlab("")+ylab(""))
      ## plot(rowMaxs(z), col= apply(z, 1, which.max)) #(z[,1]
    }
  }

  Y <- Y * gene_sd + gene_mean

  return(list("loglik" = tot, "pi" = pi, "mu" = mu, "sigma" = sigma, "beta" = beta, "lambda" = lambda, "z" = z, "Ef" = Wm, "Varf" = M, "Y" = Y, "geneM" = gene_mean, "geneSd" = gene_sd))
}
