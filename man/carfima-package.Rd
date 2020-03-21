\name{carfima-package}
\alias{carfima-package}
\docType{package}
\title{Continuous-Time Fractionally Integrated ARMA Process for Irregularly Spaced Long-Memory Time Series Data}
\description{
The R package \pkg{carfima} provides a toolbox to fit a continuous-time fractionally integrated ARMA process (CARFIMA) on univariate and irregularly spaced time series data via both frequentist and Bayesian machinery. A general-order CARFIMA(\eqn{p, H, q}) model for \eqn{p>q} is specified in Tsai and Chan (2005). It involves \eqn{p+q+2} unknown model parameters, i.e., \eqn{p} AR parameters, \eqn{q} MA parameters, Hurst parameter \eqn{H}, and process uncertainty (standard deviation) \eqn{\sigma}; see also \code{\link{carfima}}. Also, the model can account for heteroscedastic measurement errors, if the information about measurement error standard deviations is known. The package produces their maximum likelihood estimates and asymptotic uncertainties using a global optimizer called the differential evolution algorithm. It also produces posterior samples of the model parameters via Metropolis-Hastings within a Gibbs sampler equipped with adaptive Markov chain Monte Carlo. These fitting procedures, however, may produce numerical errors if \eqn{p>2}. The toolbox also contains a function to simulate discrete time series data from CARFIMA(\eqn{p, H, q}) process given the model parameters and observation times.
}
\details{
  \tabular{ll}{
    Package: \tab carfima\cr
    Type: \tab Package\cr
    Version: \tab 2.0.2\cr
    Date: \tab 2020-03-20\cr
    License: \tab GPL-2\cr
    Main functions: \tab \code{\link{carfima}}, \code{\link{carfima.loglik}}, \code{\link{carfima.sim}}\cr
  }
}
\author{
  Hyungsuk Tak, Henghsiu Tsai, and Kisung You

  Maintainer: Hyungsuk Tak <hyungsuk.tak@gmail.com>
}
\references{
    H. Tsai and K.S. Chan (2005) "Maximum Likelihood Estimation of Linear Continuous Time Long Memory Processes with Discrete Time Data," Journal of the Royal Statistical Society (Series B), 67 (5), 703-716. DOI: 10.1111/j.1467-9868.2005.00522.x

    H. Tsai and K.S. Chan (2000) "A Note on the Covariance Structure of a Continuous-time ARMA Process," Statistica Sinica, 10, 989-998.\cr Link: http://www3.stat.sinica.edu.tw/statistica/j10n3/j10n317/j10n317.htm
}
\examples{
  ##### Irregularly spaced observation time generation.

  length.time <- 30
  time.temp <- rexp(length.time, rate = 2)
  time <- rep(NA, length.time + 1)
  time[1] <- 0
  for (i in 2 : (length.time + 1)) {
    time[i] <- time[i - 1] + time.temp[i - 1]
  }
  time <- time[-1]

  ##### Data genration for CARFIMA(1, H, 0) based on the observation times. 

  parameter <- c(-0.4, 0.8, 0.2) 
  # AR parameter alpha = -0.4
  # Hurst parameter = 0.8
  # Process uncertainty (standard deviation) sigma = 0.2

  me.sd <- rep(0.05, length.time)
  # Known measurement error standard deviations 0.05 for all observations
  # If not known, remove the argument "measure.error = me.sd" in the following codes,
  # so that the default values (zero) are automatically assigned.

  y <- carfima.sim(parameter = parameter, time = time, 
                   measure.error = me.sd, ar.p = 1, ma.q = 0)  

  ##### Fitting the CARFIMA(1, H, 0) model on the simulated data for MLEs.
  \donttest{
  res <- carfima(Y = y, time = time, measure.error = me.sd, 
                 method = "mle", ar.p = 1, ma.q = 0)
  }
  # It takes a long time due to the differential evolution algorithm (global optimizer).
  # res$mle; res$se; res$AIC; res$fitted.values

  ##### Fitting the CARFIMA(1, H, 0) model on the simulated data for Bayesian inference.
  \donttest{
  res <- carfima(Y = y, time = time, measure.error = me.sd, 
                 method = "bayes", ar.p = 1, ma.q = 0, 
                 bayes.param.ini = parameter, 
                 bayes.param.scale = c(rep(0.2, length(parameter))), 
                 bayes.n.warm = 100, bayes.n.sample = 1000)
  }
  # It takes a long time because the likelihood evaluation is computationally heavy.
  # The last number of bayes.param.scale is to update sigma2 (not sigma) on a log scale.
  # hist(res$param[, 1]); res$accept; res$AIC; res$fitted.values

  ##### Computing the log likelihood of the CARFIMA(1, H, 0) model given the parameters.
  loglik <- carfima.loglik(Y = y, time = time, ar.p = 1, ma.q = 0,
                           measure.error = me.sd,
                           parameter = parameter, fitted = FALSE)
}

\keyword{package}