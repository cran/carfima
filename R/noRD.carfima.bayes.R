#' @keywords internal
#' @noRd
carfima.bayes <- function(Y, time, param.ini, param.scale, ar.p, ma.q, n.warm, n.sample){
  n.total <- n.warm + n.sample
  dimen <- length(param.ini)
  p <- ar.p
  q <- ma.q
  
  ## now main iteration runs on CPP
  tmpout = cpp_carfima_bayes(carfima.loglik.internal, Y, time, 
                             param.ini, param.scale, as.integer(p), as.integer(q), 
                             as.integer(n.warm), as.integer(n.sample))
  param.out  = tmpout$pre_param
  accept.out = tmpout$pre_accept
  
  ## return after removing n.warm samples
  out <- list(param = param.out[-c(1 : n.warm), ], 
              accept = colMeans(accept.out[-c(1 : n.warm), ]))
  return(out)
}