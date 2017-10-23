\name{carfima}			
\alias{carfima}
\title{Fitting a CARFIMA(p, H, q) model via frequentist or Bayesian machinery}
\description{
A general-order CARFIMA(\eqn{p, H, q}) model for \eqn{p>q} is
\deqn{Y_t^{(p)} -\alpha_p Y_t^{(p-1)} -\cdots- \alpha_1 Y_t = \sigma(B_{t, H}^{(1)}+\beta_1B_{t, H}^{(2)}+\cdots+\beta_q B_{t, H}^{(q+1)}),}
where \eqn{B_{t, H} = B_t^H} is the standard fractional Brownian motion, \eqn{H} is the Hurst parameter, and the superscript \eqn{(j)} indicates \eqn{j}-fold differentiation with respect to \eqn{t}; see Equation (1) of Tsai and Chan (2005) for details. The model has \eqn{p+q+2} unknown model parameters; \eqn{p} \eqn{\alpha_j}'s, \eqn{q} \eqn{\beta_j}'s, \eqn{H}, and \eqn{\sigma}.\cr

The function \code{carfima} fits the model, producing either their maximum likelihood estimates (MLEs) with their asymptotic uncertainties or their posterior samples according to its argument, \code{method}. The MLEs except \eqn{\sigma} are obtained from a profile likelihood by a global optimizer called the differential evolution algorithm on restricted ranges, i.e., (-0.99, -0.01) for each \eqn{\alpha_j}, (-10, 10) for each \eqn{\beta_j}, and (0.51, 0.99) for \eqn{H}; the MLE of \eqn{\sigma} is then deterministically computed. The corresponding asymptotic uncertainties are based on a numerical hessian matrix calculation at the MLEs (see function \code{hessian} in \pkg{numDeriv}). It also computes the Akaike Information Criterion (AIC) that is \eqn{-2}(log likelihood \eqn{-p-q-2}). The function \code{carfima} becomes numerically unstable if \eqn{p>2}, and thus it may produce numerical errors. (The original Fortran code used in Tsai and Chan (2005) does not have this numerical issue because they use a different global optimizer called \code{DBCONF} from the IMSL Fortran library.)\cr

The Bayesian approach uses independent prior distributions for the unknown model parameters; a Uniform(-0.99, -0.01) prior for each \eqn{\alpha_j}, a Normal(0, \eqn{10^4}) prior for each \eqn{\beta_j}, a Uniform(0.51, 0.99) prior for \eqn{H} for long memory process, and finally an inverse-Gamma(shape = 2.01, scale = \eqn{10^3}) prior for \eqn{\sigma^2}. Posterior propriety holds because the prior distributions are jointly proper. It also adopts appropriate proposal density functions; a truncated Normal(current state, proposal scale) between -0.99 and -0.01 for each \eqn{\alpha_j}, a Normal(current state, proposal scale) for each \eqn{\beta_j}, a truncated Normal(current state, proposal scale) between 0.51 and 0.99 for \eqn{H}, and fianlly a Normal(log(current state), proposal scale on a log scale) for \eqn{\sigma^2}, i.e., \eqn{\sigma^2} is updated on a log scale. We sample the full posterior using Metropolis within Gibbs sampler. It also adopts adaptive Markov chain Monte Carlo (MCMC) that updates the proposal scales every 100 iterations; if the acceptance rate of the most recent 100 proposals (at the end of the \eqn{i}th 100 iterationsis) smaller than 0.3 then it multiplies \eqn{\exp(-\min(0.01, 1/\sqrt{i}))} by the current proposal scale; if it is larger than 0.3 then it multiplies \eqn{\exp(\min(0.01, 1/\sqrt{i}))} by the current proposal scale. The Markov chains equipped with this adaptive MCMC converge to the stationary distribution because the adjustment factors, \eqn{\exp(\pm\min(0.01, 1/\sqrt{i}))},  approach unity as \eqn{i} goes to infinity, satisfying the diminishing adaptation condition. The function \code{carfima} becomes numerically unstable if \eqn{p>2}, and thus it may produce numerical errors. The output returns the AIC for which we evaluate the log likelihood at the posterior medians of the unknown model parameters.\cr

If the MLE-based method produces numerical errors, we recommend running the Bayesian method that is numerically more stable (though computationally more expensive).
}

\usage{carfima(Y, time, ar.p, ma.q, method = "mle", bayes.param.ini, bayes.param.scale,
        bayes.n.warm, bayes.n.sample)}

\arguments{
  \item{Y}{
    A vector of length \eqn{k} for the observed data.
  }
  \item{time}{
    A vector of length \eqn{k} for the observation times.
  }
  \item{ar.p}{
    A positive integer for the order of the AR model. \code{ar.p} must be greater than \code{ma.q}. If \code{ar.p} is greater than 2, numerical errors may occur.
  }
  \item{ma.q}{
    A non-negative integer for the order of the MA model. \code{ma.q} must be smaller than \code{ar.p}.
  }
  \item{method}{
    Either "mle" or "bayes". Method "mle" conducts the MLE-based inference, producing MLEs and asymptotic uncertainties of the model parameters. Method "bayes" draws posterior samples of the model parameters. 
  }
  \item{bayes.param.ini}{
    Only if \code{method} is "bayes". A vector of length \eqn{p+q+2} for the initial values of \eqn{p} \eqn{\alpha_j}'s, \eqn{q} \eqn{\beta_j}'s, \eqn{H}, and \eqn{\sigma} to implement a Markov chain Monte Carlo method (Metropolis within Gibbs sampler). When a CARFIMA(2, \eqn{H}, 1) model is fitted, for example, users should set five initial values of \eqn{\alpha_1},  \eqn{\alpha_2}, \eqn{\beta_1}, \eqn{H}, and \eqn{\sigma}. 
  }
  \item{bayes.param.scale}{
    Only if \code{method} is "bayes". A vector of length \eqn{p+q+2} for jumping (proposal) scales of the Metropolis steps. Each number determines how far a Metropolis step reaches out in each parameter space. Note that the last number of this vector is the jumping scale to update \eqn{\sigma^2} on a log scale. The adaptive MCMC automatically adjusts these jumping scales during the run.
  }
  \item{bayes.n.warm}{
    Only if \code{method} is "bayes". A scalar for the number of burn-ins, i.e., the number of the first few iterations to be discarded to remove the effect of initial values.
  }
  \item{bayes.n.sample}{
    Only if \code{method} is "bayes". A scalar for the number of posterior samples for each parameter.
  }
}

\details{
  The function \code{carfiam} produces MLEs, their asymptotic uncertainties, and AIC if \code{method} is "mle". It produces the posterior samples of the model parameters, acceptance rates, and AIC if \code{method} is "bayes".
}

\value{
  The outcome of \code{carfima} is composed of:
  \item{mle}{If \code{method} is "mle". Maximum likelihood estimates of the model parameters, \eqn{p} \eqn{\alpha_j}'s, \eqn{q} \eqn{\beta_j}'s, \eqn{H}, and \eqn{\sigma}. }
  \item{se}{If \code{method} is "mle". Asymptotic uncertainties (standard errors) of the MLEs.}
  \item{param}{If \code{method} is "bayes". An \eqn{m} by \eqn{(p+q+2)} matrix where \eqn{m} is the number of posterior draws (\code{bayes.n.sample}) and the columns correspond to parameters, \eqn{p} \eqn{\alpha_j}'s, \eqn{q} \eqn{\beta_j}'s, \eqn{H}, and \eqn{\sigma}.}
  \item{accept}{If \code{method} is "bayes". A vector of length \eqn{p+q+2} for the acceptance rates of the Metropolis steps.}
  \item{AIC}{For both methods. Akaike Information Criterion, -2(log likelihood \eqn{-p-q-2}). The log likelihood is evaluated at the MLEs if \code{method} is "mle" and at the posterior medians of the unknown model parameters if \code{method} is "bayes".}
  \item{fitted.values}{For both methods. A vector of length \eqn{k} for the values of \eqn{E(Y_{t_i}\vert Y_{<t_i})}, \eqn{i=1, 2, \ldots, k}, where \eqn{k} is the number of observations and \eqn{Y_{<t_i}} represents all data observed before \eqn{t_i}. Note that \eqn{E(Y_{t_1}\vert Y_{<t_1})=0}. MLEs of the model parameters are used to compute \eqn{E(Y_{t_i}\vert Y_{<t_i})}'s if \code{method} is "mle" and posterior medians of the model parameters are used if \code{method} is "bayes".}
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

  ##### Fitting the CARFIMA(1, H, 0) model on the simulated data for MLEs.
  \donttest{
  res <- carfima(Y = y, time = time, method = "mle", ar.p = 1, ma.q = 0)
  }
  # It takes a long time due to the differential evolution algorithm (global optimizer).
  # res$mle; res$se; res$AIC; res$fitted.values

  ##### Fitting the CARFIMA(1, H, 0) model on the simulated data for Bayesian inference.
  \donttest{
  res <- carfima(Y = y, time = time, method = "bayes",
                 ar.p = 1, ma.q = 0, 
                 bayes.param.ini = parameter, 
                 bayes.param.scale = c(rep(0.2, length(parameter))), 
                 bayes.n.warm = 100, bayes.n.sample = 1000)
  }
  # It takes a long time because the likelihood evaluation is computationally heavy.
  # The last number of bayes.param.scale is to update sigma2 (not sigma) on a log scale.
  # hist(res$param[, 1]); res$accept; res$AIC; res$fitted.values

  ##### Computing the log likelihood of the CARFIMA(1, H, 0) model given the parameters.

  loglik <- carfima.loglik(Y = y, time = time, ar.p = 1, ma.q = 0,
                           parameter = parameter, fitted = FALSE)
}


\keyword{methods}
