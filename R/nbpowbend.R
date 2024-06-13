#' @title  Negative-binomial-power-bend distribution functions
#'
#' @description  The PMF, CDF and RNG functions of power-bend distributed communities with negative-binomial sampling error.
#'
#' @param x,q the observed species abundances; in 16S profiling, this is the observed 16S read count
#' @param n the number of random numbers to be generated
#' @param gcn the 16S copy number of OTUs; it should be a numeric vector of the same length as `x`
#' @param eta the sampling effort measured as the expected number of reads per individual
#' @param odisp the over-dispersion parameter `r` or `size` in the negative-binomial distribution that represents sampling bias
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
#' @details `dnbpowbend` and `dtrunc_dnbpowbend` calculates the full and truncated PMF for negative-binomial-power-bend distribution
#' @details `pnbpowbend` calculates the CDF of negative-binomial-power-bend distribution. To achieve the best accuracy for SAD fitting, it calculates the upper tail of the CDF by default
#' @details `rnbpowbend` generates random number under the negative-binomial-power-bend distribution by generating random number under the power-bend distribution and then sampling under negative-binomial distribution
#' @details `constant` denotes the name of a set of values used by Gauss-Kronrod quadrature. This packages provides two sets of these values: `G25_K51`, the default; and `G7_K15`, the same as used in the GK method in `pracma::integral`
#' @details `dcon_nbpowbend` and `pconv_nbpowbend` are the integrand functions used internally by `dnbpowbend` and `pnbpowbend`
#' @details `fit_nbpowbend` fits the negative-binomial-power-bend distribution to observed SAD with maximum likelihood and a default truncation point of 0 using Nelder-Mead optimum searching algorithm. Correction for 16S rRNA gene copy number can be done by supplying the corresponding GCN to the `gcn` parameter.
#' @details The quantile function `qnbpowbend` is not implemented and thus `rnbpowbend` only generates random values under the full distribution, making it not suitable for RAD prediction.
#' @details To predict RAD, use `representative_nbpowbend` for approximated quantiles.
#' @return `dnbpowbend`,`pnbpowbend` and `dtrunc_dnbpowbend` returns a numeric vector of the (log) probability
#' @return `rnbpowbend` returns a numeric vector of the random numbers
#' @return `fit_nbpowbend` returns a `mle2` class object containing the fitted model
#' @rdname nbpowbend
dtrunc_dnbpowbend = function(x,gcn=rep(1,length(x)),trunc=0,
                             eta=1e-3,odisp,s,lambda=0.01,minlogLambda=-log(lambda),
                             exact=1e5,constant = "G25_K51",log=FALSE,debug=FALSE){
  # argument validation check and initialization
  if(eta<=0) stop("Negative-binomial sampling effort should be positive!")
  if (!missing(lambda) && !missing(minlogLambda))
    stop("specify 'lambda' or 'minlogLambda' but not both")
  lambda = exp(-minlogLambda)
  if(min(x)<=trunc) stop("Truncation should be smaller than the minimum abundance!")
  # calculate probability mass not observed
  unique.gcn = unique(gcn)
  gcn.idx = as.numeric(factor(gcn,levels = unique.gcn))
  unique.lower.sum = pnbpowbend(q=rep(trunc,length(unique.gcn)),gcn=unique.gcn,
                                eta=eta,odisp=odisp,s=s,lambda=lambda,
                                exact=exact,constant=constant,log.p=TRUE,lower.tail=FALSE)
  lower.sum = unique.lower.sum[gcn.idx]
  # initialize by the unobserved sum
  pp = -lower.sum
  # calculate truncated probability with exceptions when x==trunc+1
  if(any(x>(trunc+1))){
    these.x = x>(trunc+1)
    pp[these.x] = pp[these.x]+
      dnbpowbend(x = x[these.x],gcn=gcn[these.x],eta=eta,odisp=odisp,s=s,lambda=lambda,
                 exact=exact,constant=constant,log=TRUE)
  }
  if(any(x==(trunc+1))){
    these.x = x==(trunc+1)
    pp[these.x] = log_one_minus_x(pp[these.x]+
        pnbpowbend(q = x[these.x],gcn=gcn[these.x],eta=eta,odisp=odisp,s=s,lambda=lambda,
                   exact=exact,constant=constant,lower.tail = FALSE,log.p=TRUE))
  }
  # return truncated probability mass
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname nbpowbend
dnbpowbend = function(x,gcn=rep(1,length(x)),
                      eta=1e-3,odisp,s,lambda=0.01,minlogLambda=-log(lambda),
                      exact=1e5,constant="G25_K51",log=FALSE,debug=FALSE){
  # argument validation check and initialization
  if(eta<=0) stop("Sampling effort should be positive!")
  if (!missing(lambda) && !missing(minlogLambda))
    stop("specify 'lambda' or 'minlogLambda' but not both")
  lambda = exp(-minlogLambda)
  # degenerate to Poisson-powerbend if overdispersion of NB is small enough
  if(odisp>1e6) return(dpoipowbend(x = x,gcn = gcn,eta = eta,s = s,lambda = lambda,
                                   log = log))
  # calculation of probability mass with regard for each x and gcn pair
  pp = sapply(1:length(x),function(i){
    # localize read number xx and GCN-adjusted sampling effort etac
    xx = x[i]
    etac = eta*gcn[i]
    # calculate the log-transformed coefficient (common terms) for each
    log.coef = log(etac)*xx+log(odisp)*odisp+
      dnbinom(x = xx,size = odisp,prob = 0.5,log = TRUE)+log(2)*(odisp+xx)
    # calculate the exact probability mass given true abundance
    log.exact.pp = log_sum_x(
      dconv_nbpowbend(x=1:exact,q = xx,
                      s=s,lambda=lambda,eta=etac,odisp=odisp,log=TRUE))-
      log(etac)*(odisp+xx)
    # calculate the gap between integration approximation and series summation
    log.remainder.pp = dconv_nbpowbend(x=exact+1,q = xx,
                                       s=s,lambda=lambda,eta=etac,odisp=odisp,log=TRUE)-
      log(etac)*(odisp+xx)-log(2)
    # calculate the discriminant for the root of the first derivative of the integrand
    delta.dis = (lambda*odisp+etac*s+etac*odisp)^2-4*lambda*etac*(s-xx)*odisp
    # and find the peak of the integrand after 'exact'
    if(delta.dis>=0){
      turn.point = (sqrt(delta.dis)-(lambda*odisp+etac*s+etac*odisp))/(2*lambda*etac)
    }else{
      turn.point = 0
    }
    if(turn.point<=(exact+1)) turn.point = exact+2
    # approximate the tail part through a preset paradigm
    log.approx.pp = integral_refined_a2inf(
      fun = dconv_nbpowbend,params = list(q=xx,s=s,lambda = 1,eta = etac/lambda,odisp = odisp),
      xmin = (exact+1)*lambda,xbreak = turn.point*lambda,fold = 100,constant = constant,
      log = TRUE,silent = TRUE)-log(etac)*(odisp+xx)+log(lambda)*(s+odisp-1)
    # return the sum of each part of the series
    log.total.pp = log_sum_x(c(log.exact.pp,log.approx.pp,log.remainder.pp))+log.coef
    return(log.total.pp)
  })
  # normalize the probability mass by polylogarithm
  pp = pp-polylog_approximation(lambda = lambda,s = s,log=TRUE)
  # return the probability
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname nbpowbend
dconv_nbpowbend = function(x,q,s,lambda,eta,odisp,adjust=0,log=FALSE){
  pp = -lambda*x-(s-q)*log(x)-log(odisp+eta*x)*(odisp+q)+log(eta)*(odisp+q)-adjust
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname nbpowbend
pnbpowbend = function(q,gcn=rep(1,length(q)),
                      eta=1e-3,odisp,s,lambda=0.01,minlogLambda=-log(lambda),
                      exact=1e5,constant="G25_K51",log.p=FALSE,lower.tail=FALSE,debug=FALSE){
  # argument validation check and initialization
  if(eta<=0) stop("Sampling effort should be positive!")
  if (!missing(lambda) && !missing(minlogLambda))
    stop("specify 'lambda' or 'minlogLambda' but not both")
  lambda = exp(-minlogLambda)
  # degenerate to Poisson-powerbend if overdispersion of NB is small (odisp is large) enough
  if(odisp>1e6) return(ppoipowbend(q = q,gcn = gcn,eta = eta,s = s,lambda = lambda,
                                   exact = exact,constant = constant,lower.tail = lower.tail,log.p = log.p))
  # calculation CDF as P(q>x)
  pp = sapply(1:length(q),function(i){
    etac = eta*gcn[i]
    qq = q[i]
    ppp.exact = log_sum_x(pconv_nbpowbend(x = 1:exact,q=qq,s = s,lambda = lambda,eta = etac,odisp = odisp,log = TRUE))
    ppp.remainder = pconv_nbpowbend(x = 1+exact,q=qq,s = s,lambda = lambda,eta = etac,odisp = odisp,log = TRUE)-log(2)
    # calculate the discriminant for the root of the first derivative of the integrand
    delta.dis = (lambda*odisp+etac*s+etac*odisp)^2-4*lambda*etac*(s-qq-1)*odisp
    if(delta.dis>=0){
      turn.point = (sqrt(delta.dis)-(lambda*odisp+etac*s+etac*odisp))/(2*lambda*etac)
    }else{
      turn.point = 0
    }
    if(turn.point<=(exact+1)) turn.point = exact+2
    ppp.approx = integral_refined_a2inf(
      fun=pconv_nbpowbend,params = list(s = s, q = qq, lambda = 1, eta = etac/lambda,odisp = odisp),
      xmin = (1+exact)*lambda, xbreak = turn.point*lambda,fold = 100, constant = constant,
      silent = TRUE,log = TRUE)+log(lambda)*(s-1)
    ppp = log_sum_x(c(ppp.exact,ppp.remainder,ppp.approx))
    return(ppp)
  })-polylog_approximation(s=s,lambda=lambda,exact=exact,log=TRUE)
  # adjust lower vs upper tail
  if(lower.tail) pp = log_one_minus_x(pp)
  # return the probability
  if(log.p) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname nbpowbend
pconv_nbpowbend = function(x,q,s,lambda,eta,odisp,adjust=0,log=FALSE){
  pp = pnbinom(q = q,size = odisp,mu = x*eta,log.p = TRUE,lower.tail = FALSE)
  if(any(is.nan(pp))) pp[is.nan(pp)] = 0
  pp = pp-lambda*x-s*log(x)-adjust
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname nbpowbend
rnbpowbend = function(n,gcn=rep(1,n),eta=1e-3,odisp,s,lambda=0.01,minlogLambda=-log(lambda)){
  if (!missing(lambda) && !missing(minlogLambda))
    stop("specify 'lambda' or 'minlogLambda' but not both")
  lambda = exp(-minlogLambda)
  n0 = rpowbend(n = n,s = s,omega=lambda)
  x = rnbinom(n = n,size = odisp,mu = n0*eta*gcn)
  return(x)
}

#' @title  Fit SAD with negative-binomial-power-bend distribution
#'
#' @description  Fit a power-bend SAD compounded with negative-binomial distribution
#'
#' @param x the observed species abundances
#' @param gcn the 16S copy number of the
#' @param tol the maximum relative error of logarithm;if exceeded, a warning message will be returned
#' @details This function is cloned frminlogLambda the internal function frminlogLambda `sads` package
#' @return returns a boolean vector
#' @export
#' @rdname fit_nbpowbend
fit_nbpowbend = function(x, gcn=rep(1,length(x)),trunc=0,
                         exact=1e5,constant = "G25_K51",
                         start.value,min.odisp=1e-3,min.eta=1e-16,...){
  # argument validation check and initialization
  dots = list(...)
  if (any(x < 0) | any(!is.wholenumber(x)))
    stop("All x must be non-negative integers")
  if (!is.null(trunc)) {
    if (min(x) <= trunc)
      stop("truncation point should be lower than the lowest data value")
  }
  if (missing(start.value)){
    start.value = list(s=1,minlogLambda=6.9,eta=0.01,odisp=1)
  }
  # reorder starting values
  start.value = start.value[c("s","minlogLambda","eta","odisp")]
  # overwrite fixed starting values
  if ("fixed" %in% names(dots)){
    start.value[names(dots$fixed)] = dots$fixed
  }
  # transform starting values to apply inherent constraints
  start.value$s = log(start.value$s)
  start.value$eta = log(start.value$eta)
  start.value$odisp = log(start.value$odisp)
  # transform fixed values as well
  if ("fixed" %in% names(dots)){
    dots$fixed = start.value[names(dots$fixed)]
  }
  # setting default optimum searching method
  if (!"method" %in% names(dots)) dots$method = "Nelder-Mead"
  #reduce computation load by calculating unique abundance/gcn pairs only
  ctdf = aggregate(1:length(x),by=data.frame(abundance=x,gcn=gcn),FUN=length)
  #cat(sprintf("Fitting SAD on %d unique abundance/copy number pairs\n",dim(ctdf)[1]))
  # localize likelihood function
  if (is.null(trunc)){
    LL = function(s, minlogLambda,eta,odisp){
      pps = dnbpowbend(x = ctdf$abundance, gcn = ctdf$gcn, s = exp(s),
                       minlogLambda = minlogLambda, eta = exp(eta)+min.eta,
                       odisp = exp(odisp)+min.odisp, exact = exact,
                       constant = constant, log = TRUE)*ctdf$x
      pp = -sum(pps)
      if(is.infinite(pp)) pp = .Machine$double.xmax
      return(pp)
    }
  }
  else{
    LL = function(s, minlogLambda,eta,odisp){
      pps = dtrunc_dnbpowbend(x = ctdf$abundance, gcn = ctdf$gcn, trunc = trunc,
                              s = exp(s), minlogLambda = minlogLambda, eta = exp(eta)+min.eta,
                              odisp = exp(odisp)+min.odisp, exact = exact,
                              constant = constant, log = TRUE)*ctdf$x
      pp = -sum(pps)
      if(is.infinite(pp)) pp = .Machine$double.xmax
      return(pp)
    }
  }
  # start optimization
  result = do.call("mle2", c(list(LL, start = start.value,
                                  data = list(x = x)), dots))
  # transform parameters in optimizer back to linear scale
  result@fullcoef["eta"] = exp(result@fullcoef["eta"])+min.eta
  result@fullcoef["odisp"] = exp(result@fullcoef["odisp"])+min.odisp
  result@fullcoef["s"] = exp(result@fullcoef["s"])
  # return optimization results
  return(result)
}

