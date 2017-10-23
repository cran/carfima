\name{carfima.loglik}			
\alias{carfima.loglik}
\title{Computing the log likelihood function of a CARFIMA(p, H, q) model}
\description{
This function evaluates the log likelihood of a CARFIMA(p, H, q) model as specified in Tsai and Chan (2005). 
}

\usage{carfima.loglik(Y, time, ar.p, ma.q, parameter, fitted = FALSE)}

\arguments{
  \item{Y}{
    A vector for the \eqn{k} observed data.
  }
  \item{time}{
    A vector for the \eqn{k} observation times.
  }
  \item{ar.p}{
    A positive integer for the order of the AR model. \code{ar.p} must be greater than \code{ma.q}. If \code{ar.p} is greater than 2, numerical errors may occur for both methods.
  }
  \item{ma.q}{
    A non-negative integer for the order of the MA model. \code{ma.q} must be smaller than \code{ar.p}.
  }
  \item{parameter}{
    The values of the unknown parameters at which the log likelihood is evaluated. For example, users need to specify five values of \eqn{\alpha_1},  \eqn{\alpha_2}, \eqn{\beta_1}, \eqn{H}, and \eqn{\sigma} for CARFIMA(2, H, 1). 
  }
  \item{fitted}{
    If "TRUE", fitted values and AIC are returned. If "FALSE", a value of the log likelihood is returned. Default is "FALSE".
  }
}

\details{
  The function \code{carfiam.loglik} computes the log likelihood of a CARFIMA(p, H, q) model via the innovation algorithm whose computational cost increases linearly as the size of the data increases; see Tsai and Chan (2005) for details.
}

\value{
  The outcome of \code{carfima} is the value of the log likelihood if "fitted = FALSE" and both AIC and fitted values if "fitted = TRUE".
}

\references{
    H. Tsai and K.S. Chan (2005) "Maximum Likelihood Estimation of Linear Continuous Time Long Memory Processes with Discrete Time Data," Journal of the Royal Statistical Society (Series B), 67 (5), 703-716. DOI: 10.1111/j.1467-9868.2005.00522.x

    H. Tsai and K.S. Chan (2000) "A Note on the Covariance Structure of a Continuous-time ARMA Process," Statistica Sinica, 10, 989-998. \cr Link: http://www3.stat.sinica.edu.tw/statistica/j10n3/j10n317/j10n317.htm
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
  # process uncertainty (standard deviation) sigma = 0.2
  y <- carfima.sim(parameter = parameter, time = time, ar.p = 1, ma.q = 0)  

  ##### Computing the log likelihood of the CARFIMA(1, H, 0) model given the parameters.

  loglik <- carfima.loglik(Y = y, time = time, ar.p = 1, ma.q = 0,
                           parameter = parameter, fitted = FALSE)
}


\keyword{methods}
