#' @title  Arithmetic under log-scale
#'
#' @description  Add or subtract extremely large or small numbers in log-scale to avoid numerical issues
#'
#' @param x values to be added or subtracted
#' @param y values to be added
#' @details `log_x_plus_y` calculates z so that exp(z) = exp(x)+exp(y)
#' @details `log_sum_x` calculates z so that exp(z) = sum(exp(x))
#' @details `log_one_minus_x` calculates z so that exp(z) = 1-exp(x)
#' @details Only `log_one_minus_x` is vectorized
#' @return `log_x_plus_y` and `log_sum_x` return a single numeric value; `log_one_minus_x` returns a numeric vector with the same length as input `x`
#' @export
#' @rdname log_sum_x
log_x_plus_y = function(x,y){
  alpha = max(x,y)
  if(any(is.infinite(c(x,y)))) return(alpha)
  beta = abs(x-y)
  return(alpha+log1p(exp(-beta)))
}

#' @export
#' @rdname log_sum_x
log_sum_x = function(x){
  max.x = max(x)
  if(is.infinite(max.x)) return(max.x)
  x = x-max.x
  sum.x = sum(exp(x))
  return(log(sum.x)+max.x)
}

#' @export
#' @rdname log_sum_x
log_one_minus_x = function(x){
  res = log(1-exp(x))
  exception1 = which((abs(x)>sqrt(.Machine$double.eps))&
              (x < (log(.Machine$double.eps)/2)))
  if(length(exception1)>0){
    # Use the approximation 1-e^x = e^res = 1+res for e^x->0 or x->-Inf
    res[exception1] = -exp(x[exception1])
  }
  exception2 = which(abs(x)<=sqrt(.Machine$double.eps))
  if(length(exception2)>0){
    # Use approximation e^x = 1+x+x^2/2 for e^x -> 1
    res[exception2] = log(-x[exception2]-x[exception2]^2/2)
  }
  return(res)
}
