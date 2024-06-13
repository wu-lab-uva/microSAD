#' @title  Calculate polylogarithm for log-series and power-bend distribution
#'
#' @description  Calculate polylogarithm for log-series and power-bend with better support for base close to 1 and relevant arguments
#'
#' @param lambda the base of the numerator; `exp(-lambda)` is identical to `z` in `pracma::polylog`
#' @param s the exponent of the denominator; `s` is identical to `n` in `pracma::polylog`
#' @param exact the number of terms that is added exactly before approximation applies
#' @param log if `TRUE`, return the log-transformed value; defaults to `FALSE`
#' @details `polylog_approximation` is the wrapper function for polylogarithm function
#' @details `polylog_tail_approx` calculates the polylogarithm with better support for `z=exp(-lambda)` close to 1 compared to pracma::polylog. The tail part of the infinite sum is approximated by the continuous incomplete gamma function. Half of the `exact+1`th element in the series is supplemented to compensate the continuous approximation.
#' @return `polylog_tail_approx` returns a single numeric value
#' @export
#' @rdname polylog_approximation
#'
polylog_approximation = function(method="_tail_approx",...){
  if(method=="_approximation") method="_tail_approx"
  dots = list(...)
  return(do.call(paste0("polylog",method),dots))
}

#' @export
#' @rdname polylog_approximation
#'
polylog_tail_approx = function(lambda,s,exact = 1e5,debug = FALSE,log=FALSE){
  if(lambda<0|is.na(lambda)) stop("lambda should non-negative for convergence")
  # When s==0, polylogarithm has close-form solution z/(1-z)
  if(s==0){
    if(log) return(-lambda-log_one_minus_x(-lambda))
    return(exp(-lambda-log_one_minus_x(-lambda)))
  }
  # When s==1, polylogarithm has close-form solution -log(1-z)
  if(s==1){
    if(log) return(log(-log_one_minus_x(-lambda)))
    return(-log_one_minus_x(-lambda))
  }
  # General solution that uses incomplete gamma function as approximation
  Vexact = log_sum_x(powbend_term(x = 1:exact,s = s,lambda = lambda,log = TRUE))
  Vgap = powbend_term(x=exact+1,s=s,lambda=lambda,log=TRUE)-log(2)
  this.incgam = get_incgam(a = 1-s,x = (exact+1)*lambda,log = TRUE)
  Vtail = this.incgam + (s-1)*log(lambda)
  liz = log_sum_x(c(Vexact,Vgap,Vtail))
  if(debug) print(c(LiZ = liz,Vtail=Vtail,Vgap=Vgap,Vexact=Vexact))
  if(log) return(liz)
  return(exp(liz))
}
