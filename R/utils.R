#' @title  Integer check
#'
#' @description  Integer check for discrete distributions and other discrete functions
#'
#' @param x the numbers to be checked
#' @param tol the tolerance to distinguish real number from an integer; defaults to the square-root of the double-float precision as in the `sads` package.
#' @details The function `is.wholenumber` is cloned from the internal function in `sads` package. It is exported so that user's customized SAD functions may use it.
#' @details This function is vectorized.
#' @return `is.wholenumber` returns a Boolean vector
#' @export
#' @rdname integercheck
is.wholenumber = function(x, tol = .Machine$double.eps^0.5){
  return(abs(x - round(x)) < tol)
}

#' @title  manually calculate AICc
#'
#' @description  calculate AICc from summary statistics
#'
#' @param aic Akaike Information Criterion calculated from AIC() or other functions, can be a vector
#' @param n the number of data points, can be a vector
#' @param k the number of parameters, can be a vector
#' @return the AICc values in a numeric vector
#' @export
#' @rdname AICc_manual
AICc_manual = function(aic,n,k){
  aicc = aic+2*(k^2+k)/(n-k-1)
  return(aicc)
}

#' @title  calculate CI from Hessian matrix
#'
#' @description  Calculate CI from Hessian matrix of a minus-log-likelihood function
#'
#' @param params fitted parameters returned by optimizer
#' @param hessian the hessian matrix returned by optimizer
#' @param vcov the variance-covariance matrix returned by optimizer
#' @param alpha desired confidence level
#' @return `calculate_CI_from_Hessian_minlog` and `calculate_CI_from_variance_covariance` return a data frame listing estimated parameters and their upper and lower CI
#' @export
#' @rdname calculate_CI_from_Hessian_minlog
calculate_CI_from_Hessian_minlog = function(params,hessian,alpha=0.05){
  if(is.null(hessian)) stop("Hessian matrix is required to calculate CI.")
  inv.hessian = solve(hessian)
  params.se = sqrt(diag(inv.hessian))
  qse = qnorm(1-alpha/2)
  uCI = params + qse*params.se
  lCI = params - qse*params.se
  return(data.frame(est=params,upper=uCI,lower=lCI))
}

#' @export
#' @rdname calculate_CI_from_Hessian_minlog
calculate_CI_from_variance_covariance = function(params,vcov,alpha=0.05){
  if(is.null(vcov)) stop("Variance-covariance matrix is required to calculate CI.")
  params.se = sqrt(diag(vcov))
  qse = qnorm(1-alpha/2)
  uCI = params + qse*params.se
  lCI = params - qse*params.se
  return(data.frame(est=params,upper=uCI,lower=lCI))
}
