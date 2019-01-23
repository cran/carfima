#' Fitting a CARFIMA(p, H, q) model via frequentist or Bayesian machinery
#' 
#' A general-order CARFIMA(\eqn{p, H, q}) model for \eqn{p>q} is 
#' \deqn{Y_t^{(p)} -\alpha_p Y_t^{(p-1)} -\cdots- \alpha_1 Y_t = \sigma(B_{t, H}^{(1)}+\beta_1B_{t, H}^{(2)}+\cdots+\beta_q B_{t, H}^{(q+1)}),}
#' where \eqn{B_{t, H} = B_t^H} is the standard fractional Brownian motion, \eqn{H} is the Hurst parameter, and the 
#' superscript \eqn{(j)} indicates \eqn{j}-fold differentiation with respect to \eqn{t}; see Equation (1) of Tsai and Chan (2005) 
#' for details. The model has \eqn{p+q+2} unknown model parameters; \eqn{p} \eqn{\alpha_j}'s, \eqn{q} \eqn{\beta_j}'s, \eqn{H}, and \eqn{\sigma}.\cr
#' 
#' The function \code{carfima} fits the model, producing either their maximum likelihood estimates (MLEs) with their asymptotic uncertainties 
#' or their posterior samples according to its argument, \code{method}. The MLEs except \eqn{\sigma} are obtained from a profile likelihood 
#' by a global optimizer called the differential evolution algorithm on restricted ranges, i.e., (-0.99, -0.01) for each \eqn{\alpha_j}, 
#' (-10, 10) for each \eqn{\beta_j}, and (0.51, 0.99) for \eqn{H}; the MLE of \eqn{\sigma} is then deterministically computed. 
#' The corresponding asymptotic uncertainties are based on a numerical hessian matrix calculation at the MLEs (see function \code{hessian} 
#' in \pkg{numDeriv}). It also computes the Akaike Information Criterion (AIC) that is \eqn{-2}(log likelihood \eqn{-p-q-2}). 
#' The function \code{carfima} becomes numerically unstable if \eqn{p>2}, and thus it may produce numerical errors. 
#' (The original Fortran code used in Tsai and Chan (2005) does not have this numerical issue because they use a different global 
#' optimizer called \code{DBCONF} from the IMSL Fortran library.)\cr
#' 
#' The Bayesian approach uses independent prior distributions for the unknown model parameters; a Uniform(-0.99, -0.01) 
#' prior for each \eqn{\alpha_j}, a Normal(0, \eqn{10^4}) prior for each \eqn{\beta_j}, a Uniform(0.51, 0.99) prior for \eqn{H} 
#' for long memory process, and finally an inverse-Gamma(shape = 2.01, scale = \eqn{10^3}) prior for \eqn{\sigma^2}. 
#' Posterior propriety holds because the prior distributions are jointly proper. It also adopts appropriate proposal density functions; 
#' a truncated Normal(current state, proposal scale) between -0.99 and -0.01 for each \eqn{\alpha_j}, a Normal(current state, proposal scale) 
#' for each \eqn{\beta_j}, a truncated Normal(current state, proposal scale) between 0.51 and 0.99 for \eqn{H}, 
#' and fianlly a Normal(log(current state), proposal scale on a log scale) for \eqn{\sigma^2}, i.e., \eqn{\sigma^2} is updated 
#' on a log scale. We sample the full posterior using Metropolis within Gibbs sampler. It also adopts adaptive Markov chain Monte Carlo (MCMC) 
#' that updates the proposal scales every 100 iterations; if the acceptance rate of the most recent 100 proposals (at the end of the 
#' \eqn{i}th 100 iterationsis) smaller than 0.3 then it multiplies \eqn{\exp(-\min(0.01, 1/\sqrt{i}))} by the current proposal scale; 
#' if it is larger than 0.3 then it multiplies \eqn{\exp(\min(0.01, 1/\sqrt{i}))} by the current proposal scale. The Markov chains 
#' equipped with this adaptive MCMC converge to the stationary distribution because the adjustment factors, \eqn{\exp(\pm\min(0.01, 1/\sqrt{i}))},  
#' approach unity as \eqn{i} goes to infinity, satisfying the diminishing adaptation condition. The function \code{carfima} becomes 
#' numerically unstable if \eqn{p>2}, and thus it may produce numerical errors. The output returns the AIC for which we evaluate 
#' the log likelihood at the posterior medians of the unknown model parameters.\cr
#' 
#' If the MLE-based method produces numerical errors, we recommend running the Bayesian method that is numerically more stable (though computationally more expensive).
#' 
#' 
#' @param Y A vector of length \eqn{k} for the observed data.
#' @param time A vector of length \eqn{k} for the observation times.
#' @param ar.p A positive integer for the order of the AR model. \code{ar.p} must be greater than \code{ma.q}. If \code{ar.p} is greater than 2, numerical errors may occur.
#' @param ma.q A non-negative integer for the order of the MA model. \code{ma.q} must be smaller than \code{ar.p}.
#' @param method Either "mle" or "bayes". Method "mle" conducts the MLE-based inference, producing MLEs and asymptotic uncertainties of the model parameters. Method "bayes" draws posterior samples of the model parameters.
#' @param bayes.param.ini Only if \code{method} is "bayes". A vector of length \eqn{p+q+2} for the initial values of \eqn{p} \eqn{\alpha_j}'s, \eqn{q} \eqn{\beta_j}'s, \eqn{H}, and \eqn{\sigma} to implement a Markov chain Monte Carlo method (Metropolis within Gibbs sampler). When a CARFIMA(2, \eqn{H}, 1) model is fitted, for example, 
#' users should set five initial values of \eqn{\alpha_1},  \eqn{\alpha_2}, \eqn{\beta_1}, \eqn{H}, and \eqn{\sigma}.
#' @param bayes.param.scale Only if \code{method} is "bayes". A vector of length \eqn{p+q+2} for jumping (proposal) scales of the Metropolis steps. Each number determines how far a Metropolis step reaches out in each parameter space. Note that the last number of this vector is the jumping scale to update \eqn{\sigma^2} on a log scale. 
#' The adaptive MCMC automatically adjusts these jumping scales during the run.
#' @param bayes.n.warm Only if \code{method} is "bayes". A scalar for the number of burn-ins, i.e., the number of the first few iterations to be discarded to remove the effect of initial values.
#' @param bayes.n.sample Only if \code{method} is "bayes". A scalar for the number of posterior samples for each parameter.
#' 
#' @section Details:
#' The function \code{carfima} produces MLEs, their asymptotic uncertainties, and AIC if \code{method} is "mle". It produces the posterior samples of the model parameters, acceptance rates, and AIC if \code{method} is "bayes".
#' 
#' @return A name list composed of:
#' \describe{
#' \item{mle}{If \code{method} is "mle". Maximum likelihood estimates of the model parameters, \eqn{p} \eqn{\alpha_j}'s, \eqn{q} \eqn{\beta_j}'s, \eqn{H}, and \eqn{\sigma}. }
#' \item{se}{If \code{method} is "mle". Asymptotic uncertainties (standard errors) of the MLEs.}
#' \item{param}{If \code{method} is "bayes". An \eqn{m} by \eqn{(p+q+2)} matrix where \eqn{m} is the number of posterior draws (\code{bayes.n.sample}) and the columns correspond to parameters, \eqn{p} \eqn{\alpha_j}'s, \eqn{q} \eqn{\beta_j}'s, \eqn{H}, and \eqn{\sigma}.}
#' \item{accept}{If \code{method} is "bayes". A vector of length \eqn{p+q+2} for the acceptance rates of the Metropolis steps.}
#' \item{AIC}{For both methods. Akaike Information Criterion, -2(log likelihood \eqn{-p-q-2}). The log likelihood is evaluated at the MLEs if \code{method} is "mle" and at the posterior medians of the unknown model parameters if \code{method} is "bayes".}
#' \item{fitted.values}{For both methods. A vector of length \eqn{k} for the values of \eqn{E(Y_{t_i}\vert Y_{<t_i})}, \eqn{i=1, 2, \ldots, k}, where \eqn{k} is the number of observations and \eqn{Y_{<t_i}} represents all data observed before \eqn{t_i}. Note that \eqn{E(Y_{t_1}\vert Y_{<t_1})=0}. MLEs of the model parameters are used to compute \eqn{E(Y_{t_i}\vert Y_{<t_i})}'s if \code{method} is "mle" and posterior medians of the model parameters are used if \code{method} is "bayes".}
#' }
#' 
#' @examples 
#' \donttest{
#' ##### Irregularly spaced observation time generation.
#' 
#' length.time <- 100
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
#' #### Estimation 1 : MLE
#' 
#' res1 <- carfima(Y = y, time = time, method = "mle", ar.p = 1, ma.q = 0)
#' 
#' #### Estimation 2 : Bayes
#' 
#' res2 <- carfima(Y = y, time = time, method = "bayes", ar.p = 1, ma.q = 0, 
#' bayes.param.ini = parameter, bayes.param.scale = c(rep(0.2, length(parameter))),
#' bayes.n.warm = 100, bayes.n.sample = 1000)
#' }
#' 
#' @references 
#' \insertRef{tsai_note_2000}{carfima} 
#' 
#' \insertRef{tsai_maximum_2005}{carfima}
#' 
#' @export
carfima <- function(Y, time, ar.p, ma.q, method = c("mle","bayes"), 
                    bayes.param.ini, bayes.param.scale, 
                    bayes.n.warm, bayes.n.sample) {
  if ((ar.p>1)||(ma.q>0)){
    return(old.carfima(Y, time, ar.p, ma.q, method = "mle",
                       bayes.param.ini, bayes.param.scale,
                       bayes.n.warm, bayes.n.sample))
  } else {
    p <- ar.p
    q <- ma.q
    
    method = match.arg(method)
    if (method == "mle") {
      
      opt.res <- DEoptim::DEoptim(fn = carfima.loglik.nosigma,  
                                  control = DEoptim::DEoptim.control(itermax = 15, 
                                                                     reltol = 1e-4, steptol = 1),
                                  lower = c(rep(-0.99, ar.p), rep(-10, ma.q), 0.51), 
                                  upper = c(rep(-0.01, ar.p), rep(10, ma.q), 0.99), 
                                  Y = Y, time = time, ar.p = p, ma.q = q, sigma.hat = FALSE)
      
      #     opt.res <- DEoptim(fn = carfima.loglik.nosigma,  
      #                        control = DEoptim.control(itermax = 15, 
      #                                                  reltol = 1e-4, steptol = 1),
      #                        lower = c(rep(-0.99, ar.p), rep(-10, ma.q), 0.51), 
      #                        upper = c(rep(-0.01, ar.p), rep(10, ma.q), 0.99), 
      #                        ar.p = ar.p, ma.q = ma.q, Y = Y, time = time, 
      #                        sigma.hat = FALSE)
      
      
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
      asymp.var <- diag(-solve(numDeriv::hessian(func = carfima.loglik, 
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
    return(output)    
  }
}








# old functions in R ------------------------------------------------------
#' @keywords internal
#' @noRd
# g_1(h_(i + 1), d_j) in Appendix 3 of Tsai and Chan (Tech. report, 2005)
old.g1 <- function(h, d, H) {

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

  return(g1.res)

  # row indicates times (from h_0 to h_r, length of r + 1)
  # column indicates eigenvalues

}

# g_2(h_i, d_j) in Appendix 3 of Tsai and Chan (Tech. report, 2005)
#' @keywords internal
#' @noRd
old.g2 <- function(h, d, H) {

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

  return(as.matrix(g2.res[-(leng.h0 + 1), ]))
  # row indicates times (from h_0 to h_r)
  # column indicates eigenvalues

}

#' @keywords internal
#' @noRd
old.B.mat <- function(p, alpha) {

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

#' @keywords internal
#' @noRd
old.V.mat.sigma <- function(p, alpha, sigma) {

  B <- old.B.mat(p, alpha)
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

#' @keywords internal
#' @noRd
old.V.mat <- function(p, alpha) {

  B <- old.B.mat(p, alpha)
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

  return(V)

}

#' @keywords internal
#' @noRd
old.Gamma.Y.sigma <- function(time.lag.cov, p, A, H, beta, delta.p, sigma) {

  Gamma.Y.out <- rep(NA, length(time.lag.cov) + 1)
  C_H <- H * (2 * H - 1)
  alpha <- A[p, ]
  eig.res <- eigen(A)
  eig.val <- eig.res$values
  eig.vec <- eig.res$vectors
  eig.vec.inv <- solve(eig.vec)
  time.lag.cov.0 <- c(0, time.lag.cov)

  V.ast <- old.V.mat.sigma(p, alpha, sigma)

  g1.save <- old.g1(h = time.lag.cov, d = eig.val, H = H)
  g2.save <- old.g2(h = time.lag.cov, d = eig.val, H = H)

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

#' @keywords internal
#' @noRd
old.Gamma.Y <- function(time.lag.cov, p, A, H, beta, delta.p) {

  Gamma.Y.out <- rep(NA, length(time.lag.cov) + 1)
  C_H <- H * (2 * H - 1)
  alpha <- A[p, ]
  eig.res <- eigen(A)
  eig.val <- eig.res$values
  eig.vec <- eig.res$vectors
  eig.vec.inv <- solve(eig.vec)
  time.lag.cov.0 <- c(0, time.lag.cov)

  V.ast <- old.V.mat(p, alpha)

  g1.save <- old.g1(h = time.lag.cov, d = eig.val, H = H)
  g2.save <- old.g2(h = time.lag.cov, d = eig.val, H = H)

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

#' @keywords internal
#' @noRd
old.carfima.loglik <- function(Y, time, ar.p, ma.q, parameter, fitted = FALSE) {

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

  Gamma.Y <- Re(old.Gamma.Y.sigma(time.lag.cov = time.lag.cov, p = p, A = A, H = H,
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

#' @keywords internal
#' @noRd
old.carfima.loglik.nosigma <- function(Y, time, ar.p, ma.q, parameter, sigma.hat = FALSE) {

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

  GammaY <- Re(old.Gamma.Y(time.lag.cov = time.lag.cov, p = p, A = A, H = H,
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

#' @keywords internal
#' @noRd
old.carfima.bayes <- function(Y, time, param.ini, param.scale, ar.p, ma.q, n.warm, n.sample) {

  n.total <- n.warm + n.sample
  dimen <- length(param.ini)
  p <- ar.p
  q <- ma.q

  accept.out <- matrix(0, nrow = n.total, ncol = dimen)
  param.out <- matrix(NA, nrow = n.total, ncol = dimen)
  param.t <- param.ini

  prev.log.den <- old.carfima.loglik(parameter = param.t,
                                 Y, time, ar.p = p, ma.q = q)

  for (i in 1 : n.total) {

    for (j in 1 : p) {
      param.p <- truncnorm::rtruncnorm(1, a = -0.99, b = -0.01,
                            mean = param.t[j], sd = param.scale[j])
      param.temp <- param.t
      param.temp[j] <- param.p
      param.p.vec <- param.temp

      prop.log.den <- old.carfima.loglik(parameter = param.p.vec,
                                     Y, time, ar.p = p, ma.q = q)

      l.metro <- prop.log.den - prev.log.den
      l.hastings <- log(truncnorm::dtruncnorm(param.t[j], a = -0.99, b = -0.01,
                                   mean = param.p, sd = param.scale[j])) -
        log(truncnorm::dtruncnorm(param.p, a = -0.99, b = -0.01,
                       mean = param.t[j], sd = param.scale[j]))

      if (l.metro + l.hastings > -stats::rexp(1)) {
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
        param.p <- stats::rnorm(1, mean = param.t[p + j], sd = param.scale[p + j])
        param.temp <- param.t
        param.temp[p + j] <- param.p
        param.p.vec <- param.temp

        prop.log.den <- old.carfima.loglik(parameter = param.p.vec,
                                       Y, time, ar.p = p, ma.q = q)
        l.metro <- prop.log.den - prev.log.den +
          dnorm(param.p, mean = 0, sd = 1e4, log = TRUE) -
          dnorm(param.t[p + j], mean = 0, sd = 1e4, log = TRUE)


        if (l.metro > -stats::rexp(1)) {
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
    param.p <- truncnorm::rtruncnorm(1, a = 0.51, b = 0.99,
                          mean = param.t[p + q + 1], sd = param.scale[p + q + 1])
    param.temp <- param.t
    param.temp[p + q + 1] <- param.p
    param.p.vec <- param.temp

    prop.log.den <- old.carfima.loglik(parameter = param.p.vec,
                                   Y, time, ar.p = p, ma.q = q)
    l.metro <- prop.log.den - prev.log.den

    l.hastings <- log(truncnorm::dtruncnorm(param.t[p + q + 1], a = 0.51, b = 0.99,
                                 mean = param.p, sd = param.scale[p + q + 1])) -
      log(truncnorm::dtruncnorm(param.p , a = 0.51, b = 0.99,
                     mean = param.t[p + q + 1], sd = param.scale[p + q + 1]))

    if (l.metro + l.hastings > -stats::rexp(1)) {
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
    param.p <- exp(stats::rnorm(1, mean = 2 * log(param.t[p + q + 2]), sd = param.scale[p + q + 2]))
    param.temp <- param.t
    param.temp[p + q + 2] <- sqrt(param.p)
    param.p.vec <- param.temp

    prop.log.den <- old.carfima.loglik(parameter = param.p.vec,
                                   Y, time, ar.p = p, ma.q = q)
    l.metro <- prop.log.den +
      invgamma::dinvgamma(param.p, shape = 2.01, scale = 1e3, log = TRUE) -
      prev.log.den -
      invgamma::dinvgamma(param.t[p + q + 2]^2, shape = 2.01, scale = 1e3, log = TRUE)

    l.hastings <- log(param.p) - 2 * log(param.t[p + q + 2])

    # Accept-reject
    if (l.metro + l.hastings > -stats::rexp(1)) {
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

  }

  out <- list(param = param.out[-c(1 : n.warm), ],
              accept = colMeans(accept.out[-c(1 : n.warm), ]))

  out

}
#' @keywords internal
#' @noRd
old.carfima <- function(Y, time, ar.p, ma.q, method = "mle",
                    bayes.param.ini, bayes.param.scale,
                    bayes.n.warm, bayes.n.sample) {
  p <- ar.p
  q <- ma.q

  if (method == "mle") {

    opt.res <- DEoptim(fn = old.carfima.loglik.nosigma,
                       control = DEoptim.control(itermax = 15,
                                                 reltol = 1e-4, steptol = 1),
                       lower = c(rep(-0.99, ar.p), rep(-10, ma.q), 0.51),
                       upper = c(rep(-0.01, ar.p), rep(10, ma.q), 0.99),
                       ar.p = ar.p, ma.q = ma.q, Y = Y, time = time,
                       sigma.hat = FALSE)

    log.lik <- -opt.res$optim$bestval / 2
    MLE.wo.sigma <- as.numeric(opt.res$optim$bestmem)
    MLE.sigma <- old.carfima.loglik.nosigma(parameter = MLE.wo.sigma,
                                        Y = Y, time = time, sigma.hat = TRUE,
                                        ar.p = ar.p, ma.q = ma.q)
    temp <- old.carfima.loglik(parameter = c(MLE.wo.sigma, MLE.sigma),
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

    sample.res <- old.carfima.bayes(Y = Y, time = time,
                                param.ini = bayes.param.ini,
                                param.scale = bayes.param.scale,
                                ar.p = ar.p, ma.q = ma.q,
                                n.warm = bayes.n.warm, n.sample = bayes.n.sample)
    param.median <- apply(sample.res$param, 2, median)
    temp <- old.carfima.loglik(parameter = param.median,
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
  output
}
