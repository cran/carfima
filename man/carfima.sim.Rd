\name{carfima.sim}			
\alias{carfima.sim}
\title{Simulating a CARFIMA(p, H, q) time series}
\description{
  The funstion \code{carfima.sim} produces discrete time series data that follow a CARFIMA(p, H, q) model given the values for the model parameters and observation times.
}

\usage{carfima.sim(parameter, time, ar.p, ma.q)}

\arguments{
  \item{parameter}{
    A vector of length \eqn{p+q+2} for the generative values of the model parameters; \eqn{p} values of \eqn{\alpha_j}'s, \eqn{q} values of \eqn{\beta_j}'s, \eqn{H}, and \eqn{\sigma}.
  }
  \item{time}{
    A vector for the \eqn{k} observation times, either regularly or irregularly spaced.
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
  Hyungsuk Tak and Henghsiu Tsai
}
\examples{
  ##### Irregularly spaced observation time generation.

  length.time <- 100
  time.temp <- rexp(length.time, rate = 2)
  time <- rep(NA, length.time + 1)
  time[1] <- 0
  for (i in 2 : (length.time + 1)) {
    time[i] <- time[i - 1] + time.temp[i - 1]
  }
  time <- time[-1]

  ##### Data genration for CARFIMA(1, H, 0) based on the observation times. 

  parameter <- c(-0.4, 0.75, 0.2) 
  # AR parameter alpha = -0.4
  # Hurst parameter = 0.75
  # process uncertainty sigma = 0.2
  y <- carfima.sim(parameter = parameter, time = time, ar.p = 1, ma.q = 0)  
}


\keyword{methods}
