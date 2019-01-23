#' Computing the log likelihood function of a CARFIMA(p, H, q) model
#' 
#' @param Y A vector for the \eqn{k} observed data.
#' @param time A vector for the \eqn{k} observation times.
#' @param ar.p A positive integer for the order of the AR model. \code{ar.p} must be greater than \code{ma.q}. 
#' If \code{ar.p} is greater than 2, numerical errors may occur for both methods.
#' @param ma.q A non-negative integer for the order of the MA model. \code{ma.q} must be smaller than \code{ar.p}.
#' @param parameter The values of the unknown parameters at which the log likelihood is evaluated. 
#' For example, users need to specify five values, \eqn{\alpha_1}, \eqn{\alpha_2}, \eqn{\beta_1}, \eqn{H}, and \eqn{\sigma} 
#' for \code{CARFIMA(2,H,1)}.
#' @param fitted If \code{TRUE}, fitted values and AIC are returned as a list. \code{FALSE} flag returns a value of the log likelihood.
#' 
#' @return Either a list of fitted values(\code{fitted}) and AIC(\code{AIC}), or a numeric value of the log likelihood.
#' 
#' @section Details:
#' The function \code{carfima.loglik} computes the log likelihood of a \code{CARFIMA(p,H,q)} model via the innovation algorithm 
#' whose computational cost increases linearly as the size of the data increases. See the reference for details. 
#' 
#' @examples 
#' \donttest{
#' ##### Irregularly spaced observation time generation.
#' length.time <- 100
#' time.temp <- rexp(length.time, rate = 2)
#' time <- rep(NA, length.time + 1)
#' time[1] <- 0
#' for (i in 2 : (length.time + 1)) {
#'   time[i] <- time[i - 1] + time.temp[i - 1]
#'   }
#'   time <- time[-1]
#'
#' ##### Data genration for CARFIMA(1, H, 0) based on the observation times.
#' parameter <- c(-0.4, 0.75, 0.2)
#' # AR parameter alpha = -0.4
#' # Hurst parameter = 0.75
#' # process uncertainty (standard deviation) sigma = 0.2
#' y <- carfima.sim(parameter = parameter, time = time, ar.p = 1, ma.q = 0)
#'
#' ##### Compute
#' output = carfima::carfima.loglik(Y=y,time=time,ar.p=1,ma.q=0,parameter=parameter,fitted=TRUE)
#' }
#' 
#' @references 
#' \insertRef{tsai_note_2000}{carfima} 
#' 
#' \insertRef{tsai_maximum_2005}{carfima}
#' 
#' @export
carfima.loglik <- function(Y, time, ar.p, ma.q, parameter, fitted=FALSE){
  if ((ar.p>1)||(ma.q>0)){
    return(old.carfima.loglik(Y, time, ar.p, ma.q, parameter, fitted = FALSE))
  } else {
    tmp = cpp_carfima_loglik(Y, time, ar.p, ma.q, parameter)
    
    Gamma.Y      = tmp$Gamma_Y
    Y            = tmp$Y
    time.lag.cov = tmp$time_lag_cov
    
    
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
                                               (time[i + 1] - time[k + 1]))[1] + 1] / nu[1]
          # Gamma.Y starts with Gamma.Y(0), and then Gamma.Y(h_1), Gamma.Y(h_2), ...
          # nu starts with nu_0, and then nu_1, nu_2, ...
        } else {
          theta[i, (i - k)] <- ( Gamma.Y[which(time.lag.cov ==
                                                 (time[i + 1] - time[k + 1]))[1] + 1] -
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
      return(loglik)
    } else {
      out <- list(fitted = Y.hat, AIC = AIC)
      return(out)
    } 
  }
}


#' @keywords internal
#' @noRd
carfima.loglik.internal <- function(Y, time, ar.p, ma.q, parameter){
  tmp = cpp_carfima_loglik(Y, time, ar.p, ma.q, parameter)
  
  Gamma.Y      = tmp$Gamma_Y
  Y            = tmp$Y
  time.lag.cov = tmp$time_lag_cov
  
  
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
                                             (time[i + 1] - time[k + 1]))[1] + 1] / nu[1]
        # Gamma.Y starts with Gamma.Y(0), and then Gamma.Y(h_1), Gamma.Y(h_2), ...
        # nu starts with nu_0, and then nu_1, nu_2, ...
      } else {
        theta[i, (i - k)] <- ( Gamma.Y[which(time.lag.cov ==
                                               (time[i + 1] - time[k + 1]))[1] + 1] -
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
  
  return(loglik)
}

