#' Simulating a CARFIMA(p,H,q) time series
#' 
#' The function \code{carfima.sim} produces discrete time series data that follow a \code{CARFIMA(p,H,q)} model 
#' given the values for the model parameters and observation times.
#' 
#' @param parameter A vector of length \eqn{p+q+2} for the generative values of the model parameters; 
#' \eqn{p} values of \eqn{\alpha_j}'s, \eqn{q} values of \eqn{\beta_j}'s, \eqn{H} and \eqn{\sigma}.
#' @param time A vector for the \eqn{k} observation times, either regularly or irregularly spaced.
#' @param ar.p A scalar for the order of the AR model.
#' @param ma.q A scalar for the order of the MA model.
#' 
#' @section Details:
#' This function produces simulated discrete time series data following a \code{CARFIMA(p,H,q)} model given the 
#' values for the model parameters and observation times. It first derives a \eqn{k}-dimensional multivariate 
#' Gaussian distribution whose mean set to a vector of zeroes, where \eqn{k} is the number of observations. 
#' The covariance matrix is then filled with \eqn{Cov(Y_{t_i}, Y_{t_j})} and its closed-form formula is 
#' specified in Theorem 1(b) and 1(c) of Tsai and Chan (2005).
#' 
#' @return The outcome of \code{carfima.sim} is a vector for \eqn{k} simulated data following a \code{CARFIMA(p,H,q)} 
#' model given the values for the model parameters and observation times.
#' 
#' @examples 
#' ##### Irregularly spaced observation time generation.
#' ##### For CRAN testing, time is set to be very short.
#' 
#' length.time <- 10
#' time.temp <- rexp(length.time, rate = 2)
#' time <- rep(NA, length.time + 1)
#' time[1] <- 0
#' for (i in 2 : (length.time + 1)) {
#'   time[i] <- time[i - 1] + time.temp[i - 1]
#' }
#' time <- time[-1]
#' 
#' ##### Data genration for CARFIMA(1, H, 0) based on the observation times. 
#' 
#' parameter <- c(-0.4, 0.75, 0.2) 
#' # AR parameter alpha = -0.4
#' # Hurst parameter = 0.75
#' # process uncertainty sigma = 0.2
#' y <- carfima.sim(parameter = parameter, time = time, ar.p = 1, ma.q = 0)  
#'
#'  
#' @references 
#' \insertRef{tsai_maximum_2005}{carfima}
#' 
#' @export
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
  
  Gamma.Y <- Re(Gamma.Y.sigma.gen(time.lag.cov = time.lag.cov, p = p, A = A, H = H, 
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
  
  y.sim <- MASS::mvrnorm(1, mu = rep(0, leng.time), Sigma = cov.mat)
  return(y.sim)
}



#   -----------------------------------------------------------------------
#' @keywords internal
#' @noRd
Gamma.Y.sigma.gen <- function(time.lag.cov, p, A, H, beta, delta.p, sigma) {
  
  Gamma.Y.out <- rep(NA, length(time.lag.cov) + 1)
  C_H <- H * (2 * H - 1)
  alpha <- A[p, ]
  eig.res <- eigen(A)
  eig.val <- eig.res$values
  eig.vec <- eig.res$vectors
  eig.vec.inv <- solve(eig.vec)
  time.lag.cov.0 <- c(0, time.lag.cov)
  
  V.ast <- V.mat.sigma.gen(p, alpha, sigma)
  
  g1.save <- g1.gen(h = time.lag.cov, d = eig.val, H = H)
  g2.save <- g2.gen(h = time.lag.cov, d = eig.val, H = H)
  
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
  
  return(Gamma.Y.out)
}

#' @keywords internal
#' @noRd
V.mat.sigma.gen <- function(p, alpha, sigma) {
  
  B <- B.mat.gen(p, alpha)
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
  
  return(V)
}

#' @keywords internal
#' @noRd
g1.gen <- function(h, d, H) {
  
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
      real.part.integration <- stats::integrate(integration.real, 
                                         lower = h0[i],
                                         upper = h0[i + 1], 
                                         h.val = h0[i + 1], d = d[j], H = H)$value
      img.part.integration <- stats::integrate(integration.img, 
                                        lower = h0[i],
                                        upper = h0[i + 1], 
                                        h.val = h0[i + 1], d = d[j], H = H)$value
      
      g1.res[i + 1, j] <- exp(d[j] * (h0[i + 1] - h0[i])) * 
        g1.res[i, j] + real.part.integration + 
        sqrt(as.complex(-1)) * img.part.integration                          
    }
  }
  
  return(g1.res)
  
  # row indicates times (from h_0 to h_r, length of r + 1)
  # column indicates eigenvalues
  
}

#' @keywords internal
#' @noRd
# g_2(h_i, d_j) in Appendix 3 of Tsai and Chan (Tech. report, 2005)
g2.gen <- function(h, d, H) {
  
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
      real.part.integration <- stats::integrate(integration.real, 
                                         lower = h0.pseudo.max[i],
                                         upper = h0.pseudo.max[i + 1],
                                         h.value = h0.pseudo.max[i], d = d[j], H = H)$value
      img.part.integration <- stats::integrate(integration.img, 
                                        lower = h0.pseudo.max[i],
                                        upper = h0.pseudo.max[i + 1],
                                        h.value = h0.pseudo.max[i], d = d[j], H = H)$value
      
      g2.res[i, j] <- exp(d[j] * (h0.pseudo.max[i + 1] - h0.pseudo.max[i])) * 
        g2.res[i + 1, j] + real.part.integration + 
        sqrt(as.complex(-1)) * img.part.integration                          
    }
  }
  
  return(as.matrix(g2.res[-(leng.h0 + 1), ]))
  # row indicates times (from h_0 to h_r)
  # column indicates eigenvalues
  
}

#' @keywords internal
#' @noRd
B.mat.gen <- function(p, alpha) {
  
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
  
  return(B)
}