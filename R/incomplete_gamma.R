#' @title  Calculate incomplete gamma function
#'
#' @description  Calculate upper incomplete gamma function using recursive equation
#' @param a exponent of the power-law term
#' @param x the starting point of the infinite integration
#' @param log if `TRUE`, log-transformed value is returned; defaults to `FALSE`
#' @details `get_incgam` is an extension from `pracma::incgam` to handle situations where `a` is very large or negative
#' @details It uses the recurrence relation: incgam(x,a+1) = incgam(x,a)a + exp(-x)x^a
#' @details It is NOT vectorized for either `a` or `x`
#' @return `get_incgam` returns a single numeric value
#' @export
#' @rdname get_incgam
get_incgam = function(a,x,log=FALSE,debug=FALSE){
  # for small positive a values, call pracma::incgam directly
  if(a>=0&a<10){
    if(log) return(log(pracma::incgam(x = x,a = a)))
    return(pracma::incgam(x = x,a = a))
  }
  # for large a values iterate up from small to large
  if(a>=10){
    num.iteration = floor(a)-9
    start.a = 9+(a-floor(a))
    this.incgam = log(pracma::incgam(x = x,a = start.a))
    for(i in 1:num.iteration){
      this.incgam = log_x_plus_y(this.incgam+log(start.a),log(x)*(start.a)-x)
      start.a = start.a+1
    }
  }else{# for negative a values iterate down from the smallest non-negative value
    num.iteration = -floor(a)
    start.a = a-floor(a)
    this.incgam = log(pracma::incgam(x = x,a = start.a))
    for(i in 1:num.iteration){
      start.a = start.a-1
      powbend.term = log(x)*start.a-x
      # exception 1: when incgam yields strict zero value, output zero
      if(is.infinite(this.incgam)){
        break;
      }
      # exception 2: when incgam yields very inaccurate value, output zero
      if(this.incgam >= powbend.term){
        this.incgam = -Inf
        break;
      }
      this.incgam = powbend.term + log_one_minus_x(this.incgam - powbend.term) - log(-start.a)
    }
  }
  if(log) return(this.incgam)
  return(exp(this.incgam))
}
