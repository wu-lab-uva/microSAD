#' @title  Power-bend SAD with Poisson sampling error
#'
#' @description  The PMF, CDF, RNG and fitting function of power-bend distributed communities with Poisson sampling error.
#'
#' @param x,q the observed species abundances, non-negative integers
#' @param n the number of random numbers to be generated
#' @param gcn the 16S copy number of each species, with the same length and order of x, positive integers
#' @param eta the sampling effort, positive real number
#' @param trunc the truncation point to fit a distribution whose lower-tail is truncated; non-negative integer
#' @param s the s parameter of power-bend distribution; it is the same as parameter s in `sads` package
#' @param lambda,minglogLambda the lambda parameter of power-bend distribution; it is the same as parameter omega and oM in `sads` package
#' @param exact the number of terms to be calculated exactly before approximation applies
#' @param constant the name of Gauss-Kronrod quadrature formula; See details
#' @param log,log.p if log-transformed probability should be returned or received; default to FALSE
#' @param lower.tail if the lower-tail cumulative probabilty should be returned or received; default to FALSE for accuracy
#' @param start.value a list of named numeric values as the starting values of ML fitting
#' @param ... other parameters passed to the `optim` function
#' @param adjust log-adjustment value to avoid precision problem in integration
#' @details `dpoipowbend` and `dtrunc_dpoipowbend` calculates the full and truncated PMF for Poisson-power-bend distribution
#' @details `ppoipowbend` calculates the CDF of Poisson-power-bend distribution. To achieve the best accuracy for SAD fitting, it calculates the upper tail of the CDF by default
#' @details `rpoipowbend` generates random number under the Poisson-power-bend distribution by generating random number under the power-bend distribution and then sampling under Poisson distribution
#' @details `constant` denotes the name of a set of values used by Gauss-Kronrod quadrature. This packages provides two sets of these values: `G25_K51`, the default; and `G7_K15`, the same as used in the GK method in `pracma::integral`
#' @details `pconv_poipowbend` is the integrand function used internally by `ppoipowbend`
#' @details `fit_poipowbend` fits the Poisson-power-bend distribution to observed SAD with maximum likelihood and a default truncation point of 0 using Nelder-Mead optimum searching algorithm. Correction for 16S rRNA gene copy number can be done by supplying the corresponding GCN to the `gcn` parameter.
#' @details The quantile function `qpoipowbend` is not implemented and thus `rpoipowbend` only generates random values under the full distribution, making it not suitable for RAD prediction.
#' @details To predict RAD, use `representative_poipowbend` for approximated quantiles.
#' @return `dpoipowbend`,`ppoipowbend` and `dtrunc_dpoipowbend` returns a numeric vector of the (log) probability
#' @return `rpoipowbend` returns a numeric vector of the random numbers
#' @return `fit_poipowbend` returns a `mle2` class object containing the fitted model
#' @export
#' @rdname poipowbend
dtrunc_dpoipowbend = function(x,gcn=rep(1,length(x)),trunc=0,eta=1e-3,s,lambda=0.01,minlogLambda=-log(lambda),
                              exact=1e5,constant="G25_K51",log=FALSE,debug=FALSE){
  if(eta<0) stop("Poisson sampling effort should be non-negative!")
  if (!missing(lambda) && !missing(minlogLambda))
    stop("specify 'lambda' or 'minlogLambda' but not both")
  lambda = exp(-minlogLambda)
  if(is.nan(lambda)) stop("Parameter lambda should be non-negative!")
  if(min(x)<=trunc) stop("Truncation should be smaller than the minimum abundance!")
  if(lambda==0) return(dtrunc_poipower(x = x,gcn = gcn,trunc = trunc,eta = eta,s = s,
                                       exact = exact,log = log,debug = debug))
  # calculate probability mass not observed
  unique.gcn = unique(gcn)
  gcn.idx = as.numeric(factor(gcn,levels = unique.gcn))
  unique.lower.sum = ppoipowbend(q=rep(trunc,length(unique.gcn)),gcn=unique.gcn,
                                eta=eta,s=s,lambda=lambda,
                                exact=exact,constant=constant,log.p=TRUE,lower.tail=FALSE)
  lower.sum = unique.lower.sum[gcn.idx]
  # initialize by the unobserved sum
  pp = -lower.sum
  # calculate truncated probability with exceptions when x==trunc+1
  if(any(x>(trunc+1))){
    these.x = x>(trunc+1)
    pp[these.x] = pp[these.x]+
      dpoipowbend(x = x[these.x],gcn=gcn[these.x],eta=eta,s=s,lambda=lambda,
                 log=TRUE)
  }
  if(any(x==(trunc+1))){
    these.x = x==(trunc+1)
    pp[these.x] = log_one_minus_x(pp[these.x]+
                                    ppoipowbend(q = x[these.x],gcn=gcn[these.x],eta=eta,s=s,lambda=lambda,
                                                exact=exact,constant=constant,lower.tail = FALSE,log.p=TRUE))
  }
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname poipowbend
dpoipowbend = function(x,gcn=rep(1,length(x)),eta=1e-3,s,lambda=0.01,minlogLambda=-log(lambda),
                       log=FALSE,debug=FALSE){
  if(eta<0) stop("Poisson sampling effort should be non-negative!")
  if (!missing(lambda) && !missing(minlogLambda))
    stop("specify 'lambda' or 'minlogLambda' but not both")
  lambda = exp(-minlogLambda)
  if(is.nan(lambda)) stop("Parameter lambda should be non-negative!")
  if(lambda==0) return(dpoipower(x = x,gcn = gcn,eta = eta,s = s,
                                 log = log,debug = debug))
  pp = sapply(1:length(x),function(i){
    xx = x[i]
    if(is.infinite(xx)) return(-Inf)
    etac = gcn[i]*eta
    log.coef = dpois(x = xx,lambda = etac,log = TRUE)+etac
    log.exact.pp = polylog_approximation(lambda = lambda+etac,s = s-xx,log = TRUE)
    return(log.exact.pp+log.coef)
  })-polylog_approximation(lambda = lambda,s = s,log = TRUE)
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname poipowbend
ppoipowbend = function(q,gcn=rep(1,length(q)),
                       eta=1e-3,s,lambda=0.01,minlogLambda=-log(lambda),
                       exact=1e5,constant="G25_K51",log.p=FALSE,lower.tail=TRUE,debug=FALSE){
  if(eta<0) stop("Poisson sampling effort should be non-negative!")
  if (!missing(lambda) && !missing(minlogLambda))
    stop("specify 'lambda' or 'minlogLambda' but not both")
  lambda = exp(-minlogLambda)
  if(is.nan(lambda)) stop("Parameter lambda should be non-negative!")
  if(lambda==0) return(ppoipower(q = q,gcn = gcn,eta = eta,s = s,
                                 log = log,debug = debug))
  polylog.coef = polylog_approximation(lambda = lambda,s = s,log = TRUE)
  pp = sapply(1:length(q),function(i){
    etac = gcn[i]*eta
    qq = q[i]
    if(qq<0) return(polylog.coef)
    if(is.infinite(qq)) return(-Inf)
    ppp.exact = log_sum_x(pconv_poipowbend(x = 1:exact,q=qq,s = s,lambda = lambda,eta = etac,log = TRUE))
    ppp.remainder = pconv_poipowbend(x = 1+exact,q=qq,s = s,lambda = lambda,eta = etac,log = TRUE)-log(2)
    # calculate the discriminant for the root of the first derivative of the integrand
    delta.dis = s-qq-1
    if(delta.dis>=0){
      turn.point = 0
    }else{
      turn.point = exp(-log_one_minus_x((lambda+etac)/delta.dis))
      if(is.infinite(turn.point)) turn.point = 100*exact
    }
    if(turn.point<=(exact+1)) turn.point = exact+2
    ppp.approx = integral_refined_a2inf(
      fun=pconv_poipowbend,params = list(s = s, q = qq, lambda = 1, eta = etac/lambda),
      xmin = (1+exact)*lambda, xbreak = turn.point*lambda,fold = 100, constant = constant,
      silent = TRUE,log = TRUE)+log(lambda)*(s-1)
    ppp = log_sum_x(c(ppp.exact,ppp.remainder,ppp.approx))
    return(ppp)
  })-polylog.coef
  if(lower.tail) pp = log_one_minus_x(pp)
  if(log.p) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname poipowbend
pconv_poipowbend = function(x,q,s,lambda,eta,adjust=0,log=FALSE){
  pp = -lambda*x-log(x)*s+ppois(q = q,lambda = eta*x,lower.tail = FALSE,log.p = TRUE)-adjust
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname poipowbend
rpoipowbend = function(n,gcn=rep(1,n),eta=1e-3,s,lambda=0.01,minlogLambda=-log(lambda)){
  if (!missing(lambda) && !missing(minlogLambda))
    stop("specify 'lambda' or 'minlogLambda' but not both")
  lambda = exp(-minlogLambda)
  if(is.nan(lambda)) stop("Parameter lambda should be non-negative!")
  if(lambda==0) stop("RNG for Poisson-power has not been implemented yet!")
  n0 = r_powbend(n = n,s = s,lambda=lambda)
  x = rpois(n = n,lambda = n0*eta*gcn)
  return(x)
}

#' @export
#' @rdname poipowbend
fit_poipowbend = function(x, gcn=rep(1,length(x)),trunc=0, exact=1e5,start.value,...){
  dots = list(...)
  if (any(x < 0) | any(!is.wholenumber(x)))
    stop("All x must be non-negative integers")
  if (!is.null(trunc)) {
    if (min(x) <= trunc)
      stop("truncation point should be lower than the lowest data value")
  }
  if (missing(start.value)){
    start.value = list(s=1,minlogLambda=6.9,eta=0.01)
  }
  start.value = start.value[c("s","minlogLambda","eta")]
  if ("fixed" %in% names(dots)){
    start.value[names(dots$fixed)] = dots$fixed
  }
  start.value$s = log(start.value$s)
  start.value$eta = log(start.value$eta)
  if ("fixed" %in% names(dots)){
    dots$fixed = start.value[names(dots$fixed)]
  }
  if (!"method" %in% names(dots)) dots$method = "Nelder-Mead"
  if (dots$method == "L-BFGS-B") {
    if (!"lower" %in% names(dots))
      dots$lower = c(s = log(0.1), minlogLambda = 1,eta=log(1e-20))
    if (!"upper" %in% names(dots))
      dots$upper = c(s = log(2.999), minlogLambda = 16,eta=log(10))
  }
  ctdf = aggregate(1:length(x),by=data.frame(abundance=x,gcn=gcn),FUN=length)
  #cat(sprintf("Fitting SAD on %d unique abundance/copy number pairs\n",dim(ctdf)[1]))
  if (is.null(trunc)) {
    LL = function(s, minlogLambda,eta){
      s = exp(s)
      eta = exp(eta)
      pps = dpoipowbend(x = ctdf$abundance,gcn=ctdf$gcn,
                       s = s, minlogLambda = minlogLambda, eta = eta,log = TRUE)*ctdf$x
      pp = -sum(pps)
      if(is.infinite(pp)) pp = .Machine$double.xmax
      return(pp)
    }
  }
  else {
    LL = function(s, minlogLambda,eta){
      s = exp(s)
      eta = exp(eta)
      pps = dtrunc_dpoipowbend(x = ctdf$abundance,gcn=ctdf$gcn,
                              s = s,minlogLambda = minlogLambda, eta = eta,trunc = trunc, log = TRUE)*ctdf$x
      pp = -sum(pps)
      if(is.infinite(pp)) pp = .Machine$double.xmax
      return(pp)
    }
  }
  result = do.call("mle2", c(list(LL, start = start.value,
                                  data = list(x = x)), dots))
  result@fullcoef["eta"] = exp(result@fullcoef["eta"])
  result@fullcoef["s"] = exp(result@fullcoef["s"])
  return(result)
}
