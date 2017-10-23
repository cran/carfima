
# g_1(h_(i + 1), d_j) in Appendix 3 of Tsai and Chan (Tech. report, 2005)
g1 <- function(h, d, H) {

  h0 <- c(0, h)
  leng.h <- length(h)
  leng.d <- length(d)

  g1.res <- matrix(0, nrow = leng.h + 1, ncol = leng.d)
  # it has a row of zeros, i.e., g_1(0, d_j) = 0 for all j

  integration.real <- function(u, h.val, d, H) {
    exp(Re(d) * (h.val - u)) * cos(Im(d) * (h.val - u)) * u^(2 * H - 2)    
  }

  integration.img <- function(u, h.val, d, H) {
    exp(Re(d) * (h.val - u)) * sin(Im(d) * (h.val - u)) * u^(2 * H - 2)    
  }

  for (j in 1 : leng.d) {
    for (i in 1 : leng.h) {
      real.part.integration <- integrate(integration.real, 
                                    lower = h0[i],
                                    upper = h0[i + 1], 
                                    h.val = h0[i + 1], d = d[j], H = H)$value
      img.part.integration <- integrate(integration.img, 
                                    lower = h0[i],
                                    upper = h0[i + 1], 
                                    h.val = h0[i + 1], d = d[j], H = H)$value

      g1.res[i + 1, j] <- exp(d[j] * (h0[i + 1] - h0[i])) * 
                          g1.res[i, j] + real.part.integration + 
                          sqrt(as.complex(-1)) * img.part.integration                          
    }
  }

  g1.res

  # row indicates times (from h_0 to h_r, length of r + 1)
  # column indicates eigenvalues

}

# g_2(h_i, d_j) in Appendix 3 of Tsai and Chan (Tech. report, 2005)
g2 <- function(h, d, H) {

  leng.h <- length(h)
  h0 <- c(0, h)
  leng.h0 <- length(h0)
  h0.pseudo.max <- c(h0, h0[leng.h0] * 100)
  leng.d <- length(d)

  integration.real <- function(u, h.value, d, H) {
    exp(Re(d) * (u - h.value)) * cos(Im(d) * (u - h.value)) * u^(2 * H - 2)    
  }

  integration.img <- function(u, h.value, d, H) {
    exp(Re(d) * (u - h.value)) * sin(Im(d) * (u - h.value)) * u^(2 * H - 2)    
  }

  g2.res <- matrix(0, nrow = leng.h0 + 1, ncol = leng.d)

  for (j in 1 : leng.d) {
    for (i in leng.h0 : 1) {
      real.part.integration <- integrate(integration.real, 
                                    lower = h0.pseudo.max[i],
                                    upper = h0.pseudo.max[i + 1],
                                    h.value = h0.pseudo.max[i], d = d[j], H = H)$value
      img.part.integration <- integrate(integration.img, 
                                    lower = h0.pseudo.max[i],
                                    upper = h0.pseudo.max[i + 1],
                                    h.value = h0.pseudo.max[i], d = d[j], H = H)$value

      g2.res[i, j] <- exp(d[j] * (h0.pseudo.max[i + 1] - h0.pseudo.max[i])) * 
                      g2.res[i + 1, j] + real.part.integration + 
                          sqrt(as.complex(-1)) * img.part.integration                          
    }
  }

  as.matrix(g2.res[-(leng.h0 + 1), ])
  # row indicates times (from h_0 to h_r)
  # column indicates eigenvalues

}

B.mat <- function(p, alpha) {

  B <- matrix(0, nrow = p, ncol = p)

  for (i in 1 : p) {
    for (j in 1 : p) {
      if ((2 * j - i) <= p & (2 * j - i) >= 1) {
        B[i, j] <- (-1)^(j - i) * alpha[2 * j - i]
      } else if ((2 * j - i) == (p + 1)) {
        B[i, j] <- (-1)^(j - i - 1)
      }
    }
  }

  B  

}

V.mat.sigma <- function(p, alpha, sigma) {
  
  B <- B.mat(p, alpha)
  delta.p <- rep(0, p)
  delta.p[p] <- 1
  V.diag <- -(sigma^2 / 2) * solve(B) %*% delta.p
 
  V <- matrix(0, nrow = p, ncol = p)

  for (i in 1 : p) {
    for (j in 1 : p) {
      if ((i + j) %% 2 == 0) {
        V[i, j] <- (-1)^((i - j) / 2) * V.diag[(i + j) / 2]
      } else {
        V[i, j] <- 0
      }
    }
  }

  for (u in 1 : p) {
    V[u, u] <- V.diag[u]
  }

  V

}

V.mat <- function(p, alpha) {
  
  B <- B.mat(p, alpha)
  delta.p <- rep(0, p)
  delta.p[p] <- 1
  V.diag <- -(1 / 2) * solve(B) %*% delta.p
 
  V <- matrix(0, nrow = p, ncol = p)

  for (i in 1 : p) {
    for (j in 1 : p) {
      if ((i + j) %% 2 == 0) {
        V[i, j] <- (-1)^((i - j) / 2) * V.diag[(i + j) / 2]
      } else {
        V[i, j] <- 0
      }
    }
  }

  for (u in 1 : p) {
    V[u, u] <- V.diag[u]
  }

  V

}

Gamma.Y.sigma <- function(time.lag.cov, p, A, H, beta, delta.p, sigma) {
  
  Gamma.Y.out <- rep(NA, length(time.lag.cov) + 1)
  C_H <- H * (2 * H - 1)
  alpha <- A[p, ]
  eig.res <- eigen(A)
  eig.val <- eig.res$values
  eig.vec <- eig.res$vectors
  eig.vec.inv <- solve(eig.vec)
  time.lag.cov.0 <- c(0, time.lag.cov)

  V.ast <- V.mat.sigma(p, alpha, sigma)

  g1.save <- g1(h = time.lag.cov, d = eig.val, H = H)
  g2.save <- g2(h = time.lag.cov, d = eig.val, H = H)

  M1 <- array(NA, dim = c(p, p, length(time.lag.cov) + 1))
  for (i in 1 : (length(time.lag.cov) + 1)) {
    if (p == 1) {
      M1[, , i] <- g1.save[i, ]
    } else {
      M1[, , i] <- diag(g1.save[i, ])
    }
  }

  M2 <- array(NA, dim = c(p, p, length(time.lag.cov) + 1))
  for (i in 1 : (length(time.lag.cov) + 1)) {
    if (p == 1) {
      M2[, , i] <- g2.save[i, ]
    } else {
      M2[, , i] <- diag(g2.save[i, ])
    }
  }

  for (i in 1 : (length(time.lag.cov) + 1)) {
    Gamma.Y.out[i] <- C_H * t(beta) %*% eig.vec %*% 
    (M1[, , i] + M2[, , i] + exp(time.lag.cov.0[i] * eig.val) * M2[, , 1]) %*% 
    eig.vec.inv %*% V.ast %*% beta
  }

  Gamma.Y.out

}


Gamma.Y <- function(time.lag.cov, p, A, H, beta, delta.p) {
  
  Gamma.Y.out <- rep(NA, length(time.lag.cov) + 1)
  C_H <- H * (2 * H - 1)
  alpha <- A[p, ]
  eig.res <- eigen(A)
  eig.val <- eig.res$values
  eig.vec <- eig.res$vectors
  eig.vec.inv <- solve(eig.vec)
  time.lag.cov.0 <- c(0, time.lag.cov)

  V.ast <- V.mat(p, alpha)

  g1.save <- g1(h = time.lag.cov, d = eig.val, H = H)
  g2.save <- g2(h = time.lag.cov, d = eig.val, H = H)

  M1 <- array(NA, dim = c(p, p, length(time.lag.cov) + 1))
  for (i in 1 : (length(time.lag.cov) + 1)) {
    if (p == 1) {
      M1[, , i] <- g1.save[i, ]
    } else {
      M1[, , i] <- diag(g1.save[i, ])
    }
  }

  M2 <- array(NA, dim = c(p, p, length(time.lag.cov) + 1))
  for (i in 1 : (length(time.lag.cov) + 1)) {
    if (p == 1) {
      M2[, , i] <- g2.save[i, ]
    } else {
      M2[, , i] <- diag(g2.save[i, ])
    }
  }

  for (i in 1 : (length(time.lag.cov) + 1)) {
    Gamma.Y.out[i] <- C_H * t(beta) %*% eig.vec %*% 
    (M1[, , i] + M2[, , i] + exp(time.lag.cov.0[i] * eig.val) * M2[, , 1]) %*% 
    eig.vec.inv %*% V.ast %*% beta
  }

  Gamma.Y.out

}

carfima.loglik <- function(Y, time, ar.p, ma.q, parameter, fitted = FALSE) {

  p <- ar.p
  q <- ma.q

  if (q != 0) {
    alpha <- parameter[1 : p]
    beta <- c(1, parameter[(p + 1) : (p + q)], rep(0, p - q - 1))
    H <- parameter[p + q + 1]
    sigma <- parameter[p + q + 2]
  } else {
    alpha <- parameter[1 : p]
    beta <- c(1, rep(0, p - 1))
    H <- parameter[p + 1]
    sigma <- parameter[p + 2]
  }

  delta.p <- c(rep(0, p - 1), 1)

  Y <- Y - mean(Y)
  # making mean equal to zero in the beginning

  time.lag <- matrix(NA, nrow = length(time), ncol = length(time))

  for (j in length(time) : 2) {
    for (i in j : 2) {
      time.lag[j, i - 1] <- time[j] - time[i - 1]
    }
  }

  time.lag.cov <- sort(unique(time.lag[!is.na(time.lag)]))

  A <- matrix(0, nrow = p, ncol = p)

  if (p != 1) {
    for (i in 1 : (p - 1)) {
      A[i, i + 1] <- 1
    }
  }

  A[p, ] <- alpha

  Gamma.Y <- Re(Gamma.Y.sigma(time.lag.cov = time.lag.cov, p = p, A = A, H = H, 
                              beta = beta, delta.p = delta.p, sigma = sigma))
  # Gamma.Y(0), Gamma.Y(h_1), ..., Gamma.Y(h_r)

  nu <- rep(NA, length(time))
  Y.hat <- rep(NA, length(time))
  theta <- matrix(NA, ncol = length(time) - 1, nrow = length(time) - 1)

  nu[1] <- Gamma.Y[1]
  # nu_0 = Gamma.Y(0)
  # nu starts with nu_0, and then nu_1, nu_2, ...

  Y.hat[1] <- 0
  # Y.hat(t_1) <- 0

  for (i in 1 : (length(time) - 1)) {
    for (k in 0 : (i - 1)) {
      if (k == 0) {
        theta[i, (i - k)] <- Gamma.Y[which(time.lag.cov == 
                                     (time[i + 1] - time[k + 1])) + 1] / nu[1]
        # Gamma.Y starts with Gamma.Y(0), and then Gamma.Y(h_1), Gamma.Y(h_2), ...
        # nu starts with nu_0, and then nu_1, nu_2, ...
      } else {
        theta[i, (i - k)] <- ( Gamma.Y[which(time.lag.cov == 
                                      (time[i + 1] - time[k + 1])) + 1] -
                               sum(sapply(0 : (k - 1), function(j) { 
                                          theta[k, k - j] * theta[i, i - j] * nu[j + 1] 
                                          })) ) / nu[k + 1]
      }
    }   
    Y.hat[i + 1] <- sum(theta[i, !is.na(theta[i, ])] * (Y[i : 1] - Y.hat[i : 1]))
    nu[i + 1] <- nu[1] - sum(theta[i, !is.na(theta[i, ])]^2 * nu[i : 1])
  }

  loglik <- sum(dnorm(Y, mean = Y.hat, sd = sqrt(nu), log = TRUE))
  AIC <- -2 * (loglik - ar.p - ma.q - 2) 

  if (fitted == FALSE) {
    loglik
  } else {
    out <- list(fitted = Y.hat, AIC = AIC)
    out
  }
#  sum(log(nu)) + length(time) * log(sum((Y - Y.hat)^2 / nu))

}

carfima.loglik.nosigma <- function(Y, time, ar.p, ma.q, parameter, sigma.hat = FALSE) {

  p <- ar.p
  q <- ma.q

  if (q != 0) {
    alpha <- parameter[1 : p]
    beta <- c(1, parameter[(p + 1) : (p + q)], rep(0, p - q - 1))
    H <- parameter[p + q + 1]
#    sigma <- exp(parameter[p + q + 2])
  } else {
    alpha <- parameter[1 : p]
    beta <- c(1, rep(0, p - 1))
    H <- parameter[p + 1]
#    sigma <- exp(parameter[p + 2])
  }

  delta.p <- c(rep(0, p - 1), 1)

  Y <- Y - mean(Y)
  # making mean equal to zero at the beginning

  time.lag <- matrix(NA, nrow = length(time), ncol = length(time))

  for (j in length(time) : 2) {
    for (i in j : 2) {
      time.lag[j, i - 1] <- time[j] - time[i - 1]
    }
  }

  time.lag.cov <- sort(unique(time.lag[!is.na(time.lag)]))

  A <- matrix(0, nrow = p, ncol = p)

  if (p != 1) {
    for (i in 1 : (p - 1)) {
      A[i, i + 1] <- 1
    }
  }

  A[p, ] <- alpha

  GammaY <- Re(Gamma.Y(time.lag.cov = time.lag.cov, p = p, A = A, H = H, 
                             beta = beta, delta.p = delta.p))
  # GammaY(0), GammaY(h_1), ..., GammaY(h_r)

  nu <- rep(NA, length(time))
  Y.hat <- rep(NA, length(time))
  theta <- matrix(NA, ncol = length(time) - 1, nrow = length(time) - 1)

  nu[1] <- GammaY[1]
  # nu_0 = GammaY(0)
  # nu starts with nu_0, and then nu_1, nu_2, ...

  Y.hat[1] <- 0
  # Y.hat(t_1) <- 0

  for (i in 1 : (length(time) - 1)) {
    for (k in 0 : (i - 1)) {
      if (k == 0) {
        theta[i, (i - k)] <- GammaY[which(time.lag.cov == time[i + 1] - time[k + 1]) + 
                                     1] / nu[1]
        # GammaY starts with GammaY(0), and then GammaY(h_1), GammaY(h_2), ...
        # nu starts with nu_0, and then nu_1, nu_2, ...
      } else {
        theta[i, (i - k)] <- ( GammaY[which(time.lag.cov == time[i + 1] - time[k + 1]) +
                                       1] -
                               sum(sapply(0 : (k - 1), function(j) { 
                                          theta[k, k - j] * theta[i, i - j] * nu[j + 1] 
                                          })) ) / nu[k + 1]
      }
    }   
    Y.hat[i + 1] <- sum(theta[i, !is.na(theta[i, ])] * (Y[i : 1] - Y.hat[i : 1]))
    nu[i + 1] <- nu[1] - sum(theta[i, !is.na(theta[i, ])]^2 * nu[i : 1])
  }


#  sum(dnorm(Y, mean = Y.hat, sd = sqrt(nu), log = TRUE))
  if (sigma.hat == TRUE) {
    sigma.hat <- sqrt( mean( (Y - Y.hat)^2 / nu) )
    sigma.hat
  } else {
    sum(log(nu)) + length(time) * log(sum((Y - Y.hat)^2 / nu))
  }

}

carfima.bayes <- function(Y, time, param.ini, param.scale, ar.p, ma.q, n.warm, n.sample) {

  n.total <- n.warm + n.sample
  dimen <- length(param.ini)
  p <- ar.p
  q <- ma.q

  accept.out <- matrix(0, nrow = n.total, ncol = dimen)
  param.out <- matrix(NA, nrow = n.total, ncol = dimen)
  param.t <- param.ini

  prev.log.den <- carfima.loglik(parameter = param.t, 
                          Y, time, ar.p = p, ma.q = q)

  for (i in 1 : n.total) {

    for (j in 1 : p) {
      param.p <- rtruncnorm(1, a = -0.99, b = -0.01, 
                            mean = param.t[j], sd = param.scale[j])
      param.temp <- param.t
      param.temp[j] <- param.p
      param.p.vec <- param.temp

      prop.log.den <- carfima.loglik(parameter = param.p.vec, 
                            Y, time, ar.p = p, ma.q = q)

      l.metro <- prop.log.den - prev.log.den                    
      l.hastings <- log(dtruncnorm(param.t[j], a = -0.99, b = -0.01, 
                                   mean = param.p, sd = param.scale[j])) -
                    log(dtruncnorm(param.p, a = -0.99, b = -0.01, 
                                   mean = param.t[j], sd = param.scale[j]))
                       
      if (l.metro + l.hastings > -rexp(1)) {
        param.t <- param.p.vec
        prev.log.den <- prop.log.den
        accept.out[i, j] <- 1
      }

      if (i %% 100 == 0) {
        if(mean(accept.out[i - 99 : i, j]) > 0.3) {
          scale.adj <- exp(min(0.1, 1 / sqrt(i / 100)))
        } else if (mean(accept.out[i - 99 : i, j]) < 0.3) {
          scale.adj <- exp(-min(0.1, 1 / sqrt(i / 100)))
        } else {
          scale.adj <- 1
        }
        param.scale[j] <- param.scale[j] * scale.adj
      }

    }

    if (ma.q != 0) {

      for (j in 1 : q) {

        # ma parameter update
        param.p <- rnorm(1, mean = param.t[p + j], sd = param.scale[p + j])
        param.temp <- param.t
        param.temp[p + j] <- param.p
        param.p.vec <- param.temp

        prop.log.den <- carfima.loglik(parameter = param.p.vec, 
                                  Y, time, ar.p = p, ma.q = q)
        l.metro <- prop.log.den - prev.log.den +
                   dnorm(param.p, mean = 0, sd = 1e4, log = TRUE) -
                   dnorm(param.t[p + j], mean = 0, sd = 1e4, log = TRUE)
       
                                           
        if (l.metro > -rexp(1)) {
          param.t <- param.p.vec
          prev.log.den <- prop.log.den
          accept.out[i, p + j] <- 1
        }

        if (i %% 100 == 0) {
          if(mean(accept.out[i - 99 : i, p + j]) > 0.3) {
            scale.adj <- exp(min(0.1, 1 / sqrt(i / 100)))
          } else if (mean(accept.out[i - 99 : i, p + j]) < 0.3) {
            scale.adj <- exp(-min(0.1, 1 / sqrt(i / 100)))
          } else {
            scale.adj <- 1
          }
          param.scale[p + j] <- param.scale[p + j] * scale.adj
        }

      }
    }

    # H update
    param.p <- rtruncnorm(1, a = 0.51, b = 0.99, 
                 mean = param.t[p + q + 1], sd = param.scale[p + q + 1])
    param.temp <- param.t
    param.temp[p + q + 1] <- param.p
    param.p.vec <- param.temp

    prop.log.den <- carfima.loglik(parameter = param.p.vec, 
                            Y, time, ar.p = p, ma.q = q)
    l.metro <- prop.log.den - prev.log.den
                    
    l.hastings <- log(dtruncnorm(param.t[p + q + 1], a = 0.51, b = 0.99, 
                                 mean = param.p, sd = param.scale[p + q + 1])) - 
                  log(dtruncnorm(param.p , a = 0.51, b = 0.99, 
                      mean = param.t[p + q + 1], sd = param.scale[p + q + 1]))
                       
    if (l.metro + l.hastings > -rexp(1)) {
      param.t <- param.p.vec
      prev.log.den <- prop.log.den
      accept.out[i, p + q + 1] <- 1
    }

    if (i %% 100 == 0) {
      if(mean(accept.out[i - 99 : i, p + q + 1]) > 0.3) {
        scale.adj <- exp(min(0.1, 1 / sqrt(i / 100)))
      } else if (mean(accept.out[i - 99 : i, p + q + 1]) < 0.3) {
        scale.adj <- exp(-min(0.1, 1 / sqrt(i / 100)))
      } else {
        scale.adj <- 1
      }
      param.scale[p + q + 1] <- param.scale[p + q + 1] * scale.adj
    }

    # sigma update
    param.p <- exp(rnorm(1, mean = 2 * log(param.t[p + q + 2]), sd = param.scale[p + q + 2]))
    param.temp <- param.t
    param.temp[p + q + 2] <- sqrt(param.p)
    param.p.vec <- param.temp

    prop.log.den <- carfima.loglik(parameter = param.p.vec, 
                             Y, time, ar.p = p, ma.q = q)
    l.metro <- prop.log.den + 
               dinvgamma(param.p, shape = 2.01, scale = 1e3, log = TRUE) -
               prev.log.den -
               dinvgamma(param.t[p + q + 2]^2, shape = 2.01, scale = 1e3, log = TRUE)
                    
    l.hastings <- log(param.p) - 2 * log(param.t[p + q + 2])
                       
    # Accept-reject
    if (l.metro + l.hastings > -rexp(1)) {
      param.t <- param.p.vec
      prev.log.den <- prop.log.den
      accept.out[i, p + q + 2] <- 1
    }

    if (i %% 100 == 0) {
      if(mean(accept.out[i - 99 : i, p + q + 2]) > 0.3) {
        scale.adj <- exp(min(0.1, 1 / sqrt(i / 100)))
      } else if (mean(accept.out[i - 99 : i, p + q + 2]) < 0.3) {
        scale.adj <- exp(-min(0.1, 1 / sqrt(i / 100)))
      } else {
        scale.adj <- 1
      }
      param.scale[p + q + 2] <- param.scale[p + q + 2] * scale.adj
    }

    param.out[i, ] <- param.t

    if (i %% 100 == 0) {
      print(paste(i, "iterations done:", Sys.time()))
    }

  }

  out <- list(param = param.out[-c(1 : n.warm), ], 
              accept = colMeans(accept.out[-c(1 : n.warm), ]))

  out
  
}

carfima.sim <- function(parameter, time, ar.p, ma.q) {

  p <- ar.p
  q <- ma.q

  if (q != 0) {
    alpha <- parameter[1 : p]
    beta <- c(1, parameter[(p + 1) : (p + q)], rep(0, p - q - 1))
    H <- parameter[p + q + 1]
    sigma <- parameter[p + q + 2]
  } else {
    alpha <- parameter[1 : p]
    beta <- c(1, rep(0, p - 1))
    H <- parameter[p + 1]
    sigma <- parameter[p + 2]
  }

  delta.p <- c(rep(0, p - 1), 1)
  leng.time <- length(time)

  time.lag <- matrix(NA, nrow = leng.time, ncol = leng.time)

  for (j in leng.time : 2) {
    for (i in j : 2) {
      time.lag[j, i - 1] <- time[j] - time[i - 1]
    }
  }

  time.lag.cov <- sort(unique(time.lag[!is.na(time.lag)]))

  A <- matrix(0, nrow = p, ncol = p)

  if (p != 1) {
    for (i in 1 : (p - 1)) {
      A[i, i + 1] <- 1
    }
  }

  A[p, ] <- alpha

  Gamma.Y <- Re(Gamma.Y.sigma(time.lag.cov = time.lag.cov, p = p, A = A, H = H, 
                             beta = beta, delta.p = delta.p, sigma = sigma))
  # Gamma.Y(0), Gamma.Y(h_1), ..., Gamma.Y(h_r)

  cov.mat <- matrix(NA, nrow = leng.time, ncol = leng.time)

  diag(cov.mat) <- Gamma.Y[1]
  # Cov(Y_i, Y_i) for all i
  for (i in 1 : (leng.time - 1)) {
    for (j in (i + 1) : leng.time) {
        temp <- Gamma.Y[which(time.lag.cov == time[j] - time[i])[1] + 1]
        cov.mat[i, j] <- cov.mat[j, i] <- temp
        # Gamma.Y starts with Gamma.Y(0), and then Gamma.Y(h_1), Gamma.Y(h_2), ...
    }   
  }

  y.sim <- mvrnorm(1, mu = rep(0, leng.time), Sigma = cov.mat)
  y.sim

}


carfima <- function(Y, time, ar.p, ma.q, method = "mle", 
                    bayes.param.ini, bayes.param.scale, 
                    bayes.n.warm, bayes.n.sample) {

  print(Sys.time())

  p <- ar.p
  q <- ma.q

  if (method == "mle") {

    opt.res <- DEoptim(fn = carfima.loglik.nosigma,  
                       control = DEoptim.control(itermax = 15, 
                                 reltol = 1e-4, steptol = 1),
                       lower = c(rep(-0.99, ar.p), rep(-10, ma.q), 0.51), 
                       upper = c(rep(-0.01, ar.p), rep(10, ma.q), 0.99), 
                       ar.p = ar.p, ma.q = ma.q, Y = Y, time = time, 
                       sigma.hat = FALSE)

    log.lik <- -opt.res$optim$bestval / 2
    MLE.wo.sigma <- as.numeric(opt.res$optim$bestmem)
    MLE.sigma <- carfima.loglik.nosigma(parameter = MLE.wo.sigma, 
                               Y = Y, time = time, sigma.hat = TRUE,
                               ar.p = ar.p, ma.q = ma.q)
    temp <- carfima.loglik(parameter = c(MLE.wo.sigma, MLE.sigma), 
                                          Y = Y, time = time, 
                                          ar.p = ar.p, ma.q = ma.q, 
                                          fitted = TRUE)
    fitted.values <- temp$fitted
    AIC <- temp$AIC
    asymp.var <- diag(-solve(hessian(func = carfima.loglik, 
                                     x = c(MLE.wo.sigma, MLE.sigma),
                                     Y = Y, time = time, 
                                     ar.p = ar.p, ma.q = ma.q)))

  } else {

    sample.res <- carfima.bayes(Y = Y, time = time, 
                  param.ini = bayes.param.ini, 
                  param.scale = bayes.param.scale, 
                  ar.p = ar.p, ma.q = ma.q, 
                  n.warm = bayes.n.warm, n.sample = bayes.n.sample)
    param.median <- apply(sample.res$param, 2, median)
    temp <- carfima.loglik(parameter = param.median, 
                                          Y = Y, time = time, 
                                          ar.p = ar.p, ma.q = ma.q, 
                                          fitted = TRUE)    
    fitted.values <- temp$fitted
    AIC <- temp$AIC
  }

  if (method == "bayes") {
    output <- list(param = sample.res$param, 
                   accept = sample.res$accept, AIC = AIC,
                   fitted.values = fitted.values)
  } else {
    output <- list(mle = c(MLE.wo.sigma, MLE.sigma),
                   se = sqrt(asymp.var), AIC = AIC, 
                   fitted.values = fitted.values)
  }

  print(Sys.time())  

  output
  
}

