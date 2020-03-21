\name{carfima.sim}			
\alias{carfima.sim}
\title{Simulating a CARFIMA(p, H, q) time series}
\description{
  The funstion \code{carfima.sim} produces discrete time series data that follow a CARFIMA(p, H, q) model given the values for the model parameters and observation times.
}

\usage{carfima.sim(parameter, time, measure.error, ar.p, ma.q)}

\arguments{
  \item{parameter}{
    A vector of length \eqn{p+q+2} for the generative values of the model parameters; \eqn{p} values of \eqn{\alpha_j}'s, \eqn{q} values of \eqn{\beta_j}'s, \eqn{H}, and \eqn{\sigma}.
  }
  \item{time}{
    A vector for the \eqn{k} observation times, either regularly or irregularly spaced.
  }
  \item{measure.error}{
    (Optional) A vector for the \eqn{k} measurement error standard deviations, if such information is available (especially for astronomical applications). A vector of zeros is automatically assigned, if nothing is specified.
  }
  \item{ar.p}{
    A scalar for the order of the AR model.
  }
  \item{ma.q}{
    A scalar for the order of the MA model.
  }
}

\details{
This function produces simulated discrete time series data following a CARFIMA(\eqn{p, H, q}) model given the values for the model parameters and observation times. It first derives a \eqn{k}-dimensional multivariate Gaussian distribution whose mean set to a vector of zeros, where \eqn{k} is the number of observations. The covariance matrix is filled with Cov(\eqn{Y_{t_i}}, \eqn{Y_{t_j}}) and its closed-form formula is specified in Theorem 1(b) and 1(c) of Tsai and Chan (2005).
}

\value{
  The outcome of \code{carfima.sim} is a vector for \eqn{k} simulated data following a CARFIMA(\eqn{p, H, q}) model given the values for the model parameters and observation times.
}

\references{
    H. Tsai and K.S. Chan (2005) "Maximum Likelihood Estimation of Linear Continuous Time Long Memory Processes with Discrete Time Data," Journal of the Royal Statistical Society (Series B), 67 (5), 703-716. DOI: 10.1111/j.1467-9868.2005.00522.x
}

\author{
  Hyungsuk Tak, Henghsiu Tsai, Kisung You
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
}


\keyword{methods}
