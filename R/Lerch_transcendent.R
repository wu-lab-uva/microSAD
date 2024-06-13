#' @title  Calculate Lerch transcendent for log-series and power-bend distribution
#'
#' @description  Calculate Lerch transcendent for log-series and power-bend with better support for base close to 1
#'
#' @param lambda the base of the numerator; `exp(-lambda)` is identical to `z` in `pracma::polylog`
#' @param s the exponent of the denominator; `s` is identical to `n` in `pracma::polylog`
#' @param exact the number of terms that is added exactly before approximation applies
#' @param x the index of element in the log-series or power-bend after which the sum is calculated
#' @param log if `TRUE`, return the log-transformed value; defaults to `FALSE`
#' @details The Lerch transcendent function is particularly useful for calculating power-bend CDF and the approximated q-function 
#' @details `Lerch_transcendent_approximation` is the wrapper function for Lerch transcendent function
#' @details `Lerch_transcendent_tail_approx` approximates the Lerch transcendent. The tail part of the infinite sum of the series is approximated by the continuous incomplete gamma function. Half of the `exact+1`th element in the series is supplemented to compensate the continuous approximation.
#' @return `Lerch_transcendent_tail_approx` returns a single numeric value
#' @export
#' @rdname Lerch_transcendent_approximation
#' 
Lerch_transcendent_approximation = function(method="_tail_approx",...){
  if(method=="_approximation") method="_tail_approx"
  dots = list(...)
  return(do.call(paste0("Lerch_transcendent",method),dots))
}

#' @export
#' @rdname Lerch_transcendent_approximation
#'
Lerch_transcendent_tail_approx = function(x=1,lambda,s,exact = 1e5,log=FALSE,debug=FALSE){
  Vexact = log_sum_x(powbend_term(x = x+(0:exact),s = s,lambda = lambda,log = TRUE))
  Vgap = powbend_term(x=x+exact+1,s=s,lambda=lambda,log=TRUE)-log(2)
  Vtail = get_incgam(a = 1-s,x = (x+exact+1)*lambda,log = TRUE)+(s-1)*log(lambda)
  liz = log_sum_x(c(Vexact,Vgap,Vtail))
  if(debug) print(c(LiZ = liz,Vtail=Vtail,Vgap=Vgap,Vexact=Vexact))
  if(log) return(liz)
  return(exp(liz))
}
