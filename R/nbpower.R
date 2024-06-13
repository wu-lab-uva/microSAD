#' @title  Negative-binomial-power-law distribution functions
#'
#' @description  Negative-binomial-power-law distribution functions for SAD fitting
#'
#' @param x,q the observed species abundances, non-negative integers
#' @param gcn the 16S copy number of each species, with the same length and order of x, positive integers
#' @param eta the sampling effort, positive real number
#' @param trunc the truncation point to fit a distribution whose lower-tail is truncated; non-negative integer
#' @param s exponent of the denominator in power-law distribution; it is the same as parameter `s` in the `sads` package
#' @param exact the number of terms to be calculated exactly before approximation applies
#' @param constant the name of Gauss-Kronrod quadrature constants; see `Details`
#' @param log,log.p if log-transformed probability should be returned or supplied; default to `FALSE`
#' @param lower.tail if the lower-tail cumulative probability should be returned or received; default to `FALSE` for accuracy of large quantiles
#' @param adjust log-adjustment value to avoid precision problem in integration
#' @param start.value a list of named numeric values as the starting values of ML fitting
#' @param upper the upper limit of parameters in the ML fitting
#' @param min.odisp,min.s padding of parameter values in ML fitting
#' @param ... other parameters passed to the `mle2` function
#' @details `dnbpower` is the probability mass function; `dtrunc_dnbpower` is the truncated probability mass function;`pnbpower` is the cumulative distribution function.
#' @details `dconv_nbpower` and `pconv_nbpower` are the integrand function for `dnbpower` and `pnbpower`, respectively; despite their role as internal functions, they are exported for user-created functions.
#' @details `constant` can be either a character value naming the specific element to use in the data list `quadgk_constant`, or a user-supplied list of numeric vectors containing the node and weight values.
#' @details `fit_nbpower` is the function to fit a negative-binomial power-law distribution using ML.
#' @return `dnbpower`, `dtrunc_dnbpower`, `pnbpower` return a numeric vector the same length as `x` or `q`
#' @return `fit_nbpower` returns a `mle2` class object
#' @export
#' @rdname nbpower
#'
dtrunc_dnbpower = function(x, gcn=rep(1,length(x)), trunc=0, eta, odisp, s,
                           exact=1e5, constant="G25_K51", log=FALSE, debug=FALSE){
  if(eta < 0) stop("Poisson sampling effort should be non-negative!")
  x=round(x)
  if(min(x) <= trunc) stop("Truncation should be smaller than the minimum abundance!")
  # calculate probability mass not observed
  unique.gcn = unique(gcn)
  gcn.idx = as.numeric(factor(gcn, levels = unique.gcn))

  unique.lower.sum = pnbpower(q=rep(trunc, length(unique.gcn)),
                              eta=eta, s=s,odisp=odisp,
                              exact=exact, log=TRUE, lower.tail=FALSE)
  lower.sum = unique.lower.sum[gcn.idx]

  # calculate truncated probability with exceptions when x == trunc + 1
  pp = dnbpower(x = x, eta=eta, odisp= odisp, s=s,log=TRUE)-lower.sum
  if(any(x==(trunc+1))){
    these.x = which(x==(trunc+1))
    these.x = these.x[pp[these.x]>-log(2)]
    if(length(these.x)>0){
      pp[these.x] = log_one_minus_x(pnbpower(q = x[these.x],gcn=gcn[these.x],eta=eta,s=s, odisp=odisp,
                                             exact=exact,constant=constant,lower.tail = FALSE,log=TRUE)-
                                      lower.sum[these.x])
    }
  }
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname nbpower
#'
dconv_nbpower=function(x, q, eta, odisp, s, log=FALSE, adjust = 0, debugging = FALSE){
  pp =  dnbinom(x=q, size = odisp, mu = eta*x, log = TRUE) - s*log(x) - adjust
  if (anyNA(pp)){
    these.p = which(is.na(pp))
    pp[these.p] = -Inf
  }
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname nbpower
#'
pconv_nbpower=function(x, q, eta, odisp, s, log=FALSE, adjust = 0 , debugging= FALSE){
  pp =  pnbinom(q=q, size = odisp, mu = eta*x, lower.tail=FALSE, log = TRUE) - s*log(x) - adjust
  if (anyNA(pp)){
    these.p = which(is.na(pp))
    pp[these.p] = -Inf
  }
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname nbpower
#'
dnbpower = function(x, gcn=rep(1,length(x)), eta, odisp, s,
                    exact=1e5, constant="G25_K51",
                    log = FALSE,debug=FALSE){
  if(s<=1) return(rep(NaN,length(x)))
  # degenerate to Poisson-power if overdispersion of NB is small enough
  if(odisp>1e6) return(dpoipower(x = x,gcn = gcn,eta = eta, s = s, exact=exact, log = log))
  x = floor(x)
  zeta.coef = log(VGAM::zeta(s))
  pp = sapply(1:length(x), function(i){
    xx = x[i]
    if(xx<0) return(-Inf)
    if(is.infinite(xx)) return(-Inf)
    etac = eta*gcn[i]
    ppp.exact = log_sum_x(dconv_nbpower(1:exact,q=xx,s=s,eta = etac, odisp=odisp, log = TRUE))
    ppp.remainder = dconv_nbpower(x = 1+exact,q=xx,s=s,eta = etac, odisp=odisp, log = TRUE)-log(2)
    # calculate the discriminant for the root of the first derivative of the integrand
    delta.dis = s-xx-1
    if(delta.dis>=0){
      turn.point = 0
    }else{
      turn.point = exp(-log_one_minus_x(etac/delta.dis))
      if(is.infinite(turn.point)) turn.point = 100*exact
    }
    if(turn.point<=(exact+1)) turn.point = exact+2
    ppp.tail = integral_refined_a2inf_power(
      fun=dconv_nbpower, s = s, q = xx, eta = etac, odisp=odisp,
      xmin = 1+exact, xbreak = turn.point, fold = 100,
      cutoff = -19, cutoff.type=c(1,2),silent = TRUE,log = TRUE)
    ppp.full = log_sum_x(c(ppp.exact,ppp.remainder,ppp.tail))
    return(ppp.full)
  }) - zeta.coef
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname nbpower
#'
pnbpower = function(q, gcn=rep(1,length(q)), eta, odisp, s,
                    exact=1e5, constant="G25_K51",
                    log.p = FALSE, lower.tail=TRUE, debug=FALSE){
  if(s<=1) return(rep(NaN,length(q)))
  # degenerate to Poisson-power if overdispersion of NB is small (odisp is large) enough
  if(odisp>1e6) return(ppoipower(q = q, gcn = gcn, eta = eta, s = s, exact = exact,
                                 constant = constant, lower.tail = lower.tail, log.p = log.p))
  q = floor(q)
  zeta.coef = log(VGAM::zeta(s))
  pp = sapply(1:length(q), function(i){
    qq = q[i]
    if(qq<0) return(zeta.coef)
    if(is.infinite(qq)) return(-Inf)
    etac = eta*gcn[i]
    ppp.exact = log_sum_x(pconv_nbpower(1:exact,q=qq,eta = etac, odisp=odisp,s=s,log = TRUE))
    ppp.remainder = pconv_nbpower(x = 1+exact,q=qq,s = s,eta = etac,odisp=odisp, log = TRUE)-log(2)
    # calculate the discriminant for the root of the first derivative of the integrand
    delta.dis = s-qq-1
    if(delta.dis>=0){
      turn.point = 0
    }else{
      turn.point = exp(-log_one_minus_x(etac/delta.dis))
      if(is.infinite(turn.point)) turn.point = 100*exact
    }
    if(turn.point<=(exact+1)) turn.point = exact+2
    ppp.tail = integral_refined_a2inf_power(
      fun=pconv_nbpower, s = s, q = qq, eta = etac, odisp=odisp,
      xmin = 1+exact, xbreak = turn.point, fold = 100,
      cutoff = -19, cutoff.type=c(1,2),silent = TRUE,log = TRUE)
    ppp.full = log_sum_x(c(ppp.exact,ppp.remainder,ppp.tail))
    return(ppp.full)
  }) - zeta.coef
  if(lower.tail) pp = log_one_minus_x(pp)
  if(log.p) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname nbpower
#'
fit_nbpower = function (x, gcn=rep(1,length(x)), trunc=0, start.value, upper = c(s=20), min.odisp = 1e-4,min.s=1e-10,...)
{
  dots = list(...)
  if (any(x <= 0) | any(!is.wholenumber(x)))
    stop("All x must be positive integers")
  if (!is.null(trunc)){
    if (min(x) <= trunc)
      stop("truncation point should be lower than the lowest data value")
  }
  if (missing(start.value)) {
    gamahat = function(ga, xvec) eq = -sum(log(xvec)) -
        VGAM::zeta(ga, deriv = 1) * length(xvec)/VGAM::zeta(ga)
    shat = log(uniroot(gamahat, interval = c(1.01, upper['s']), xvec = x)$root)
    etahat =log(1)
    odisphat = log(1)
  }
  else {
    shat = log(start.value$s -1 )
    etahat = log(start.value$eta)
    odisphat = log(start.value$odisp)
  }
  ctdf = aggregate(1:length(x),by=data.frame(abundance=x, gcn= gcn),FUN=length)
  if (is.null(trunc)) {
    LL = function(s,eta, odisp) {
      s = 1+min.s+exp(s)
      eta = exp(eta)
      odisp = min.odisp+exp(odisp)
      pps = dnbpower(x=ctdf$abundance, gcn=ctdf$gcn, s=s, eta=eta, odisp= odisp, log = TRUE)*ctdf$x
      pp = -sum(pps)
      if(is.infinite(pp)) pp = .Machine$double.xmax
      return(pp)}
  } else{
    LL = function(s,eta,odisp) {
      s = 1+min.s+exp(s)
      eta = exp(eta)
      odisp = min.odisp+exp(odisp)
      pps = dtrunc_dnbpower(x=ctdf$abundance, gcn=ctdf$gcn, trunc = trunc, s=s, eta=eta, odisp= odisp, log = TRUE)*ctdf$x
      pp = -sum(pps)
      if(is.infinite(pp)) pp = .Machine$double.xmax
      return(pp)}
  }
  try({
    result = do.call("mle2", c(list(LL, start = list(s = shat, eta=etahat, odisp = odisphat),
                                    data = list(x = x),method = 'L-BFGS-B' ),
                               dots))
  })
  if (!exists('result')){
    result = do.call("mle2", c(list(LL, start = list(s = shat, eta=etahat,odisp = odisphat),
                                    data = list(x = x),method = 'Nelder-Mead'),
                               dots))
  } else if (result@details$convergence !=0) {
    cat("L-BFGS-B method didnt converge, nelder-mead tried instead\n")
    result = do.call("mle2", c(list(LL, start = list(s = shat, eta=etahat),
                                    data = list(x = x),method = 'Nelder-Mead'),
                               dots))
  }
  result@fullcoef["s"] = 1+min.s+exp(result@fullcoef["s"])
  result@fullcoef["eta"] = exp(result@fullcoef["eta"])
  result@fullcoef["odisp"] = min.odisp+exp(result@fullcoef["odisp"])
  return(result)
}
