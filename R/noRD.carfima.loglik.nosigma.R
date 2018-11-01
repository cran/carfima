#' @keywords internal
#' @noRd
carfima.loglik.nosigma <- function(Y, time, ar.p, ma.q, parameter, sigma.hat = FALSE){
  tmp = cpp_carfima_loglik_nosigma(Y, time, ar.p, ma.q, parameter)
  
  GammaY       = tmp$GammaY
  Y            = tmp$Y
  time.lag.cov = tmp$time_lag_cov
  
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
        theta[i, (i - k)] <- GammaY[which(time.lag.cov == time[i + 1] - time[k + 1])[1] + 
                                      1] / nu[1]
        # GammaY starts with GammaY(0), and then GammaY(h_1), GammaY(h_2), ...
        # nu starts with nu_0, and then nu_1, nu_2, ...
      } else {
        theta[i, (i - k)] <- ( GammaY[which(time.lag.cov == time[i + 1] - time[k + 1])[1] +
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
    return(sigma.hat)
  } else {
    return(sum(log(nu)) + length(time) * log(sum((Y - Y.hat)^2 / nu)))
  }
}