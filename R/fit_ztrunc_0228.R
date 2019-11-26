
# fit-truncated normal given sigma solve mu
LL <- function(par, xbar, Sx, delta = 0) {
  sigma <- par[1]

  mu <- xbar + (Sx - sigma^2) / (xbar - delta)

  return(log(sigma) + (Sx + (xbar - mu)^2) / 2 / sigma^2 + pnorm((mu - delta) / sigma, log.p = T))
}

# LL1 <- function(par, xbar, Sx, r01, p, delta = 0)  #
# {
#   mu = par[1]
#   sigma = par[2]
#   phi =  pnorm( (delta-mu)/sigma)
#   log(sigma) + (Sx + (xbar - mu)^2)/2/sigma^2 - r01*log(1 - p + p*phi)
#
# }

# fit zero-inflated censored normal with p fixed (at either p_min or p_max), given sigma solve mu
LLwZ <- function(par, xbar, Sx, r01, p, delta = 0) #
{
  sigma <- par[1]

  mu <- xbar + (Sx - sigma^2) / (xbar - delta)

  phi <- pnorm((delta - mu) / sigma)
  log(sigma) + (Sx + (xbar - mu)^2) / 2 / sigma^2 - r01 * log(1 - p + p * phi)
}

# p_max is the dropout rate for integrating bulk
ztruncnorm <- function(dat, cutoff = 0.1, iter = 200, s_upper = 4, s_lower = 0.1, p_min = 0.6, p_max = 1) {
  n <- ncol(dat)
  dat_pos <- dat
  dat_pos[dat_pos <= cutoff] <- NA
  mu <- rowMeans(dat_pos, na.rm = T)
  var <- rowVars(dat_pos, na.rm = T)
  n1 <- rowSums(!is.na(dat_pos))
  var <- var * (n1 - 1) / n1

  if (length(p_max) == 1) p_max <- rep(p_max, nrow(dat))

  # first fit truncated normal distribution using only positive entries
  # write the log-likelihood in terms of only std and optimize by Brent
  param0 <- foreach(g = 1:nrow(dat_pos), .combine = rbind) %dopar% {
    # if(is.na(mu[g]) | is.na(var[g])) return(c(NA,NA))
    if (n1[g] < 4) {
      return(c(NA, NA))
    } # return(c(0,cutoff/2))
    # if(mu[g] - cutoff > 3.5*sqrt(var[g])) return(c(mu[g], sqrt(var[g])))
    tryCatch(
      {
        res <- optim(sqrt(var[g]), LL, xbar = mu[g], Sx = var[g], delta = cutoff, method = "Brent", upper = s_upper, lower = s_lower, control = list(maxit = iter))
        if (res$convergence != 0) {
          return(c(NA, NA))
        } else {
          c(mu[g] + (var[g] - res$par[1]^2) / (mu[g] - cutoff), res$par[1])
        }
      },
      error = function(err) {
        c(NA, NA)
      }
    )
  }

  # refit if estimated p_g is too large or too small
  n0 <- n - n1
  r01 <- n0 / n1
  p1 <- n1 / n
  p11 <- pnorm((param0[, 1] - cutoff) / param0[, 2]) # is na when too few nonzero
  pg <- p1 / p11 ## 1 - (p0 - p01)

  ind <- which(pg > p_max | is.na(pg)) # p11 == na or p1 = p11 = 0
  pg[ind] <- p_max[ind]
  ind2 <- which(pg < p_min)
  pg[ind2] <- p_min
  param1 <- foreach(g = c(ind, ind2), .combine = rbind) %dopar% {
    if (is.na(var[g])) {
      return(c(-1, 0.5))
    } # cutoff/2
    tryCatch(
      {
        res <- optim(par = sqrt(var[g]), fn = LLwZ, xbar = mu[g], Sx = var[g], r01 = r01[g], p = pg[g], delta = cutoff, method = "Brent", upper = s_upper, lower = s_lower, control = list(maxit = iter))
        # res <- optim(par=c(mu[g], sqrt(var[g])), fn=LL1, xbar = mu[g], Sx = var[g], r01 = r01[g], p = pg[g], delta = cutoff, control = list(maxit = iter))
        if (res$convergence != 0) {
          return(c(NA, NA))
        } else {
          c(mu[g] + (var[g] - res$par[1]^2) / (mu[g] - cutoff), res$par[1])
          # res$par
        }
      },
      error = function(err) {
        c(NA, NA)
      }
    )
  }

  param0[c(ind, ind2), ] <- param1

  return(list("param" = param0, "pg" = pg))
}
