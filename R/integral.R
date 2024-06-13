#' @title  Integration with Gauss-Kronrod Quadrature
#'
#' @description  A wrapper function for integration through Gauss-Kronrod (quadgk) or tanh-sinh quadrature (quadinf)
#'
#' @param fun integrand, univariate (vectorized) function; the same as in `pramca::integral`
#' @param xmin,xmax endpoints of the integration interval; the same as in `pramca::integral`
#' @param constant a character name for the Gauss-Kronrod system or a list of user-supplied numeric vectors; see `Details` for more information
#' @param no_intervals the number of subdivisions at start; the same as in `pracma::integral`
#' @param random logical; shall the length of subdivisions be random; the same as in `pracma::integral`
#' @param tol relative tolerance
#' @param silent logical; whether to output warnings and informative messages
#' @param ... additional parameters to be passed to the integrand function; the same as in `pracma::integral`
#' @details `integral_refined` is a modified version of `pracma::integral` to allow the use of higher-order GK nodes for finite integration by the `constant` parameter.
#' @details The `constant` parameter takes either a character name from names(quadgk_constant) or a list of user-provided node constants.
#' @details If provided by the user, `constant` should be a list of four elements with the name `K_node`,`K_weight`,`G_node`,`G_weight`.
#' @return `integral_refined` returns a single numeric value for the integration.
#' @seealso [quadgk_refined()]
#' @export
#' @rdname integral_refined
integral_refined = function (fun, xmin, xmax,
                             constant = "G25_K51",no_intervals = 8, random = FALSE,
                             tol=1e-8, silent=FALSE, ...){
  # argument validation check and initialization
  stopifnot(is.numeric(xmin), length(xmin) == 1,
            is.numeric(xmax), length(xmax) == 1)
  no_intervals = max(1, floor(no_intervals))
  fun = match.fun(fun)
  f = function(x) fun(x, ...)
  if (length(f(xmin)) > 1 || length(f(xmax)) > 1) {
    stop("Function 'fun' is array-valued! Use 'quadv'.\n")
  }
  if (length(f(c(xmin, xmax))) != 2) {
    if(!silent) cat("Warning: Function 'fun' is not vectorized!\n")
    f = Vectorize(f)
  }
  if (xmin == xmax) return(0)
  if(is.character(constant)){
    constant = microSAD::quadgk_constant[[constant]]
    if(is.null(constant)){
      if(!silent) cat("Warning: Supplied name for GK nodes and weights not found, using G25_K51 instead.\n")
      constant = microSAD::quadgk_constant$G25_K51
    }
  }
  # use quadinf for infinite integration
  if (is.infinite(xmin) || is.infinite(xmax)) {
    if(!silent) cat("For infinite domains Gauss integration is applied!\n")
    Q = pracma::quadinf(f, xmin, xmax, tol = tol)$Q
    return(Q)
  }
  # subdivide interval
  if (random) {
    xs = c(xmin, (xmax - xmin) * sort(runif(no_intervals - 1)) + xmin, xmax)
  }
  else {
    xs = pracma::linspace(xmin, xmax, no_intervals + 1)
  }
  # calculate integral for each subdivision
  Q = 0
  for(i in 1:no_intervals){
    Q = Q + quadgk_refined(f, xs[i], xs[i + 1], constant = constant, tol = tol)
  }
  return(Q)
}

#' @title  Adaptive Gauss-Kronrod Quadrature with user-determined nodes and weights
#'
#' @description  A modified version of `pracma::quadgk` to use user-determined nodes and weights for Gauss-Kronrod quadrature
#'
#' @param f the same as in `pracma::quadgk`: integrand as function; needs to be vectorized.
#' @param a,b the same as in `pracma::quadgk`: endpoints of the integration interval.
#' @param constant a list of node and weight values; see `Details` for more information.
#' @param tol the same as in `pracma::quadgk`: relative tolerance.
#' @param ... the same as in `pracma::quadgk`: additional parameters to be passed to the function `f`.
#' @details `quadgk_refined` is extended from `pracma::quadgk` function to use user-determined G(n)_K(2n+1) Gauss-Kronrod quadrature.
#' @details `constant` has to be a named list of 4 elements: `$K_node`,`$K_weight`,`$G_node`,`$G_weight`
#' @details `$K_node` and `K_weight` are numeric vectors of the nodes and weights for the Kronrod rule, with a length of 2n+1
#' @details `$G_node` and `G_weight` are numeric vectors of the nodes and weights for the Gauss rule, with a length of n
#' @details `microMETER` package includes two sets of node and weight values in exported data `quadgk_constant`, in which `$G25_K51` is used as the default set for `constant`.
#' @return `quadgk_refined` returns a single numeric value for the integration.
#' @seealso [integral_refined()]
#' @export
#' @rdname quadgk_refined
quadgk_refined = function(f, a, b, constant=microSAD::quadgk_constant$G25_K51,
                          tol = .Machine$double.eps^0.5, ...){
  # argument validation check and initilizations
  stopifnot(is.numeric(a), length(a) == 1, is.numeric(b), length(b) == 1)
  eps = .Machine$double.eps
  fun = match.fun(f)
  f = function(x) fun(x, ...)
  if (a == b)
    return(0)
  else if (a > b)
    return(-1 * quadgk_refined(f, b, a, constant = constant,tol = tol))
  # localize adaptive GK quadrature function
  .gkadpt = function(f, a, b, tol = tol){
    #localize node position, weight and integrand values
    xK = 0.5 * ((b - a) * constant$K_node + b + a)
    xG = 0.5 * ((b - a) * constant$G_node + b + a)
    QG = sum(constant$G_weight * f(xG)) * (b - a)/2
    QK = sum(constant$K_weight * f(xK)) * (b - a)/2
    # estimate the error and return estimated integral within tolerance
    if (!is.finite(QG) || !is.finite(QK)) {
      warning("Infinite or NA function value encountered.")
      return(QK)
    }
    else if (abs(QK - QG) < tol) {
      return(QK)
    }
    # or terminate iteration even if tolerance is not reached
    else if (abs(b - a) < 16 * eps) {
      warning("Minimum step size reached; singularity possible.")
      return(Q2)
    }
    # further divide the interval until tolerance is reached
    Q2 = .gkadpt(f, (a + b)/2, b, tol = tol)
    Q1 = .gkadpt(f, a, (a + b)/2, tol = tol)
    return(Q1 + Q2)
  }
  # estimate the integral value and return
  return(.gkadpt(f, a, b, tol = tol))
}

#' @title  Specialized infinite integration by expanding intervals
#'
#' @description  A specialized infinite integral function from a to Inf for use in probability functions
#'
#' @param fun the function for the integrand
#' @param params the parameters for the integrand
#' @param xmin the starting point of integration
#' @param xbreak the end point of the first segmented integration
#' @param fold the coefficient of exponential expansion of the interval of integration
#' @param constant the constants used in Gauss-Kronrod quadrature
#' @param no_intervals the number of subdivisions at start; the same as in `pracma::integral`
#' @param random logical; shall the length of subdivisions be random; the same as in `pracma::integral`
#' @param tol relative tolerance
#' @param silent logical; whether to output warnings and informative messages
#' @param cutoff the logged threshold of integration termination
#' @param cutoff.type the metrics used for intergration termination, see `Details`
#' @param nCheck the number of subdivisions of a finite interval for normalization of integrand values
#' @param log logical; whether to return the logged value of the integral
#' @details `integral_refined_a2inf` is a specialized protocol of integration for distribution functions like `dnbpowbend` for desired numerical accuracy in a reasonable time.
#' @details It is designed to minimize numerical instability associated with `pracma::quadinf` for some integrands by dividing the infinite interval to exponentially expanding finite intervals.
#' @details Three metrics are used to determine whether such expansion may stop:
#' @details 1. if the integral value (finite+infinite) remains stable across iterations
#' @details 2. if the infinite integral value becomes negligible to the sum of finite integrals
#' @details 3. if the increment of finite integral value becomes negligible
#' @details Not all 3 metrics are useful (and sometimes one may be harmful) for certain integrands, and user may choose which metric to use by supplying a numeric vector of the metric indices through `cutoff.type`.
#' @return `integral_refined_a2inf` returns a single numeric value for the integration.
#' @export
#' @rdname integral_a2inf
integral_refined_a2inf = function(fun, params, xmin, xbreak, fold=10,
                                  constant = "G25_K51",no_intervals = 8, random = FALSE,tol=1e-8, silent=FALSE,
                                  cutoff=-19,cutoff.type=c(1,2,3),nCheck=100,log=FALSE){
  this.min = xmin
  this.max = xbreak
  checkers = seq(from=this.min,to=this.max,length.out=nCheck)
  this.adjust = max(do.call(fun,args = c(list(x=checkers,log=TRUE),params)))
  fin.integral = log(
    do.call(integral_refined,
            args = c(list(fun = fun,xmin = this.min,xmax = this.max,
                        constant=constant,silent = silent,
                        adjust=this.adjust,log=FALSE),params)))+this.adjust
  last.increment = fin.integral
  inf.integral = log(
    do.call(integral_refined,
            args = c(list(fun = fun,xmin = this.max,xmax = Inf,
                          constant=constant,silent = silent,
                          adjust=this.adjust,log=FALSE),params)))+this.adjust
  total.integral = log_x_plus_y(fin.integral,inf.integral)
  this.delta = 1
  while(this.delta>cutoff){
    this.min = this.max
    this.max = fold*this.max
    checkers = seq(from=this.min,to=this.max,length.out=nCheck)
    this.adjust = max(do.call(fun,args = c(list(x=checkers,log=TRUE),params)))
    fin.increment = log(do.call(integral_refined,
                                args = c(list(fun = fun,xmin = this.min,xmax = this.max,
                                              constant=constant,silent = silent,
                                              adjust=this.adjust,log=FALSE),params)))+this.adjust
    fin.integral = log_x_plus_y(fin.increment,fin.integral)
    inf.integral = log(
      do.call(integral_refined,
              args = c(list(fun = fun,xmin = this.max,xmax = Inf,
                            constant=constant,silent = silent,
                            adjust=this.adjust,log=FALSE),params)))+this.adjust
    last.integral = total.integral
    total.integral = log_x_plus_y(fin.integral,inf.integral)
    this.delta = max(c(log(abs(total.integral-last.integral)),
                       inf.integral-fin.integral,
                       fin.increment-last.increment)[cutoff.type],na.rm = TRUE)
    last.increment = fin.increment
  }
  if(log) return(total.integral)
  return(exp(total.integral))
}

#' @title  Specialized infinite integration of power-law based SAD by expanding intervals
#'
#' @description  A specialized infinite integral function from a to Inf for use in power-law based SAD probability functions
#'
#' @param fun the function for the integrand
#' @param q,s,eta,odisp the parameters of the power-law-based integrand
#' @param xmin the starting point of integration
#' @param xbreak the end point of the first segmented integration
#' @param fold the coefficient of exponential expansion of the interval of integration
#' @param constant the constants used in Gauss-Kronrod quadrature
#' @param no_intervals the number of subdivisions at start; the same as in `pracma::integral`
#' @param random logical; shall the length of subdivisions be random; the same as in `pracma::integral`
#' @param tol relative tolerance
#' @param silent logical; whether to output warnings and informative messages
#' @param cutoff the logged threshold of integration termination
#' @param cutoff.type the metrics used for intergration termination, see `Details`
#' @param nCheck the number of subdivisions of a finite interval for normalization of integrand values
#' @param log logical; whether to return the logged value of the integral
#' @details `integral_refined_a2inf_power` is a specialized protocol of integration for distribution functions like `poipower` for desired numerical accuracy in a reasonable time.
#' @details It is designed to minimize numerical instability associated with `pracma::quadinf` by dividing the infinite interval to exponentially expanding finite intervals.
#' @details Two metrics are used to determine whether such expansion may stop:
#' @details 1. if the integral value (finite+infinite) remains stable across iterations
#' @details 2. if the finite integral values comply to a geometric series so analytic solutions of their sum may be used
#' @details User may choose which metric to use by supplying a numeric vector of the metric indices through `cutoff.type`.
#' @return `integral_refined_a2inf_power` returns a single numeric value for the integration.
#' @export
#' @rdname integral_a2inf_power
integral_refined_a2inf_power = function(fun, q, s, eta, odisp, xmin, xbreak, fold=10,
                                  constant = "G25_K51",no_intervals = 8, random = FALSE,tol=1e-8, silent=FALSE,
                                  cutoff=-19,cutoff.type=c(1,2),nCheck=100,log=FALSE){
  this.min = xmin
  this.max = xbreak
  if(missing(odisp)){
    static.params = list(q=q,s=s)
  }else{
    static.params = list(q=q,s=s,odisp=odisp)
  }
  checkers = seq(from=this.min/this.max,to=1,length.out=nCheck)
  this.adjust = max(do.call(fun,args = c(static.params,list(eta=eta*this.max,x=checkers,log=TRUE))))
  fin.integral = log(
    do.call(integral_refined,
            args = c(static.params,
                     list(fun = fun, eta=eta*this.max,
                          xmin = this.min/this.max,xmax = 1,constant=constant,
                          silent = silent,adjust=this.adjust,log=FALSE))))-log(this.max)*(s-1)+this.adjust
  fin.increment = fin.integral
  last.increment = fin.integral
  inf.integral = log(
    do.call(integral_refined,
            args = c(static.params,
                     list(fun = fun, eta=eta*this.max,
                          xmin = 1,xmax = Inf,constant=constant,
                          silent = silent,adjust=this.adjust,log=FALSE))))-log(this.max)*(s-1)+this.adjust
  total.integral = log_x_plus_y(fin.integral,inf.integral)
  this.delta = 1
  last.ratio = fin.increment-last.increment
  while(this.delta>cutoff){
  
    this.min = this.max
    this.max = fold*this.max
    
    if (is.infinite(this.min)){
      if(log) return(total.integral)
      return(exp(total.integral))
    }
    checkers = seq(from=this.min/this.max,to=1,length.out=nCheck)
    this.adjust = max(do.call(fun,args = c(static.params,list(eta=eta*this.max,x=checkers,log=TRUE))))
    fin.increment = log(
      do.call(integral_refined,
              args = c(static.params,
                       list(fun = fun, eta=eta*this.max,
                            xmin = this.min/this.max,xmax = 1,constant=constant,
                            silent = silent,adjust=this.adjust,log=FALSE))))-log(this.max)*(s-1)+this.adjust
    fin.integral = log_x_plus_y(fin.increment,fin.integral)
    increment.ratio = fin.increment-last.increment
    if(increment.ratio>=0){
      inf.integral = log(
        do.call(integral_refined,
                args = c(static.params,
                         list(fun = fun, eta=eta*this.max,
                              xmin = 1,xmax = Inf,constant=constant,
                              silent = silent,adjust=this.adjust,log=FALSE))))-log(this.max)*(s-1)+this.adjust
    }else{
      inf.integral = fin.increment + increment.ratio - log_one_minus_x(increment.ratio)
    }
    last.integral = total.integral
    total.integral = log_x_plus_y(fin.integral,inf.integral)
    this.delta = max(c(log(abs(total.integral-last.integral)),
                       log(abs(increment.ratio-last.ratio)))[cutoff.type],na.rm = TRUE)
    last.increment = fin.increment
    last.ratio = increment.ratio
  }
  if(log) return(total.integral)
  return(exp(total.integral))
}

#' @title  Specialized infinite integration of lognormal based SAD by expanding intervals
#'
#' @description  A specialized infinite integral function from 0 to Inf for use in lognormal based SAD probability functions
#'
#' @param fun the function for the integrand
#' @param params the parameters for the integrand
#' @param fin.intervals the logged breaking points for the finite intervals
#' @param constant the constants used in Gauss-Kronrod quadrature
#' @param no_intervals the number of subdivisions at start; the same as in `pracma::integral`
#' @param random logical; shall the length of subdivisions be random; the same as in `pracma::integral`
#' @param tol relative tolerance
#' @param silent logical; whether to output warnings and informative messages
#' @param cutoff the logged threshold of integration termination; not used for lognormal-based integrands
#' @param nCheck the number of subdivisions of a finite interval for normalization of integrand values
#' @param log logical; whether to return the logged value of the integral
#' @details `integral_refined_a2inf_lnorm` is a specialized protocol of integration for distribution functions like `dpoilog` for desired numerical accuracy in a reasonable time.
#' @details It is designed to minimize numerical instability associated with `pracma::quadinf` by dividing the interval to smaller intervals.
#' @return `integral_refined_a2inf_lnorm` returns a single numeric value for the integration.
#' @export
#' @rdname integral_02inf_lnorm
integral_refined_02inf_lnorm = function(fun, params, fin.intervals=qnorm(seq(0.0005,0.9995,length.out=9)),
                                        constant = "G25_K51",no_intervals = 8, random = FALSE,tol=1e-8, silent=FALSE,
                                        cutoff=-19,nCheck=100,log=FALSE){
  #
  this.adjust = do.call(fun,args = c(list(x=max(fin.intervals),log=TRUE),params))
  if(is.infinite(this.adjust)){
    upper.inf.integral = -Inf
  }else{
    upper.inf.integral = log(
      do.call(integral_refined,
              args = c(list(fun = fun,xmin = max(fin.intervals),xmax = Inf,
                            constant=constant,silent = silent,
                            adjust=this.adjust,log=FALSE),params)))+this.adjust
  }
  #
  this.adjust = do.call(fun,args = c(list(x=min(fin.intervals),log=TRUE),params))
  if(is.infinite(this.adjust)){
    lower.inf.integral = -Inf
  }else{
    lower.inf.integral = log(
      do.call(integral_refined,
              args = c(list(fun = fun,xmin = -Inf,xmax = min(fin.intervals),
                            constant=constant,silent = silent,
                            adjust=this.adjust,log=FALSE),params)))+this.adjust
  }
  #
  fin.intervals = lapply(1:(length(fin.intervals)-1),function(i){
    return(list(xmin=fin.intervals[i],xmax=fin.intervals[i+1]))
  })
  #
  fin.integral = sapply(fin.intervals,function(fi){
    checkers = seq(from=fi$xmin,to=fi$xmax,length.out=nCheck)
    this.adjust = max(do.call(fun,args = c(list(x=checkers,log=TRUE),params)))
    if(is.infinite(this.adjust)) return(-Inf)
    this.integral = log(
      do.call(integral_refined,
              args = c(list(fun = fun,constant=constant,silent = silent,
                            adjust=this.adjust,log=FALSE),fi,params)))+this.adjust
    return(this.integral)
  })
  #
  total.integral = log_sum_x(c(fin.integral,lower.inf.integral,upper.inf.integral))
  if(log) return(total.integral)
  return(exp(total.integral))
}
