#' @title  Negative-binomial log-normal distribution functions
#'
#' @description  The probability mass and random number generator functions for the negative-binomial log-normal distribution
#'
#' @param x,q the observed species abundances; in 16S profiling, this is the observed 16S read count
#' @param n the number of random numbers to be generated
#' @param gcn the 16S copy number of OTUs; it should be a numeric vector of the same length as `x`
#' @param odisp the over-dispersion parameter `r` or `size` in the negative-binomial distribution that represents sampling bias
#' @param mu the log mean of the log-normal distribution
#' @param sig the log standard deviation of the log-normal distribution
#' @param trunc the truncation cutoff
#' @param constant the node and weight values for Gauss-Kronrod quadrature formula; see `Details` for more information
#' @param log,log.p logical; if `TRUE`, the log-transformed value will be returned; defaults to `FALSE`
#' @details `dnblog` and `dtrunc_dnblog` are the (truncated) probability mass function; `pnblog` is the cumulative distribution function; `rnblog` is the random number generator function.
#' @details `dconv_nblog` and `pconv_nblog` are the integrand functions for `dnblog` and `pnblog`, respectively; despite their role as internal functions, it is exported for user-created functions.
#' @details In negative-binomial log-normal distribution, the sampling effort `eta`, or the mean of the negative-binomial distribution, has perfect collinearity with the `mu` in the log-normal distribution. Therefore, `eta` is degenerated and fixed to 1 internally (before GCN correction).
#' @details `constant` can be either a character value naming the specific element to use in the data list `quadgk_constant`, or a user-supplied list of numeric vectors containing the node and weight values.
#' @details The quantile function `qnblog` is not implemented and thus `rnblog` only generates random values under the full distribution, making it not suitable for RAD prediction.
#' @details For predicting RAD, use `representative_nblog` for approximated quantiles.
#' @return `dnblog`,`dtrunc_dnblog`, `pnblog`, `dconv_nblog` and `pconv_nblog` return a numeric vector the same length as `x` or `q`
#' @return `rnblog` return a numeric vector of length `n`
#' @return `fit_nblog` return a `mle2` class object of the fitted model
#' @export
#' @rdname nblog
#'
dnblog = function(x,gcn=rep(1,length(x)),mu,sig,odisp,log=FALSE,constant="G25_K51"){
  # defer to poilog when overdispersion is small (odisp large) enough
  if(odisp>1e6){
    return(d_poilog(x=x,gcn=gcn,mu = mu, sig = sig,log = log,constant = constant))
  }
  pp = sapply(1:length(x),function(i){
    xx = x[i]
    if(xx<0) return(-Inf)
    if(is.infinite(xx)) return(-Inf)
    muc = mu + log(gcn[i])
    #
    peak.est = max(log(xx),muc)
    peak.det = function(y){(muc-(odisp+1)*sig^2)*exp(y)-y*(odisp+exp(y))+odisp*((xx-1)*sig^2+muc)}
    try({peak.est = uniroot(f = peak.det, lower = min(muc,log(xx+1))-1,
                            upper = peak.est+1)$root},silent = TRUE)
    fin.intervals = peak.est+sig*seq(-5,5,1)
    fin.intervals = fin.intervals[is.finite(fin.intervals)]
    #
    ppp = integral_refined_02inf_lnorm(
      fun = dconv_nblog,params = list(q=xx,mu=muc,sig=sig,odisp=odisp),
      fin.intervals = fin.intervals,constant=constant,log=TRUE,silent = TRUE)
    return(ppp)
  })
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname nblog
#'
dconv_nblog = function(x,q,mu,sig,odisp,log=FALSE,adjust=0){
  pp =  dnorm(x = x,mean = mu,sd = sig,log = TRUE)+
    dnbinom(x = q,size = odisp,mu = exp(x),log = TRUE)-
    adjust
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname nblog
#'
dtrunc_dnblog = function(x,gcn=rep(1,length(x)),mu,sig,odisp,trunc=0,constant="G25_K51",log=FALSE,debug=FALSE){
  if(min(x)<=trunc) stop("Truncation should be smaller than the minimum abundance!")
  unique.gcn = unique(gcn)
  gcn.idx = as.numeric(factor(gcn,levels = unique.gcn))
  unique.lower.sum = pnblog(q = rep(trunc,length(unique.gcn)),gcn=unique.gcn,
                            mu = mu, sig = sig, odisp = odisp,
                            constant = constant,log.p = TRUE,lower.tail = FALSE)
  lower.sum = unique.lower.sum[gcn.idx]
  # initialize by the unobserved sum
  pp = -lower.sum
  # calculate truncated probability with exceptions when x==trunc+1
  if(any(x>(trunc+1))){
    these.x = x>(trunc+1)
    pp[these.x] = pp[these.x]+dnblog(x = x[these.x],gcn=gcn[these.x],
                                     mu=mu,sig=sig,odisp=odisp,log=TRUE)
  }
  if(any(x==(trunc+1))){
    these.x = x==(trunc+1)
    pp[these.x] = log_one_minus_x(
      pp[these.x]+pnblog(q = x[these.x],gcn=gcn[these.x],mu=mu,sig=sig,odisp=odisp,
                         constant=constant,lower.tail = FALSE,log.p=TRUE))
  }
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname nblog
#'
pnblog = function(q,gcn=rep(1,length(q)),mu,sig,odisp,log.p=FALSE,lower.tail=FALSE,constant="G25_K51"){
  # degenerate to poilog when overdispersion is small (odisp value large) enough
  if(odisp>1e6){
    return(p_poilog(q = q, gcn = gcn, mu = mu, sig = sig,
                    log.p = log.p, lower.tail=lower.tail, constant=constant))
  }
  pp = sapply(1:length(q),function(i){
    qq = q[i]
    if(qq<0) return(0)
    if(is.infinite(qq)) return(-Inf)
    muc = mu + log(gcn[i])
    #
    peak.est = max(log(qq),muc)
    peak.det = function(y){(muc-(odisp+1)*sig^2)*exp(y)-(y)*(odisp+exp(y))+odisp*(qq*sig^2+muc)}
    try({peak.est = uniroot(f = peak.det, lower = min(muc,log(qq+1))-1,
                            upper = peak.est+1)$root},silent = TRUE)
    fin.intervals = peak.est+sig*seq(-5,5,1)
    fin.intervals = fin.intervals[is.finite(fin.intervals)]
    fin.intervals = unique(round(fin.intervals,digits = 8))
    #
    ppp = integral_refined_02inf_lnorm(
      fun = pconv_nblog,params = list(q=qq,mu=muc,sig=sig,odisp=odisp),
      fin.intervals = fin.intervals,constant=constant,log = TRUE,silent = TRUE)
    return(ppp)
  })
  if(lower.tail) pp = log_one_minus_x(pp)
  if(log.p) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname nblog
#'
pconv_nblog = function(x,q,mu,sig,odisp,log=FALSE,adjust=0){
  pp = pnbinom(q = q,size = odisp,mu = exp(x),log.p = TRUE,lower.tail = FALSE)
  if(any(is.nan(pp))) pp[is.nan(pp)] = 0
  pp = pp + dnorm(x = x,mean = mu,sd = sig,log = TRUE) - adjust
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname nblog
#'
rnblog = function(n,mu,sig,odisp){
  rnbinom(n = n,size = odisp,mu = rlnorm(n = n,meanlog = mu,sdlog = sig))
}

#' @export
#' @rdname nblog
#'
fit_nblog = function(x,gcn=rep(1,length(x)),trunc=0,start.value, constant="G25_K51",min.odisp=0.001,min.sig=1e-5,...){
  dots = list(...)
  if (any(x < 0) | any(!is.wholenumber(x)))
    stop("All x must be non-negative integers")
  if (!is.null(trunc)) {
    if (min(x) <= trunc)
      stop("truncation point should be lower than the lowest data value")
  }
  if (missing(start.value)){
    start.value = list(mu=0,sig=1,odisp=1)
  }
  start.value = start.value[c("mu","sig","odisp")]
  if ("fixed" %in% names(dots)){
    start.value[names(dots$fixed)] = dots$fixed
  }
  start.value$sig = log(start.value$sig)
  start.value$odisp = log(start.value$odisp)
  if ("fixed" %in% names(dots)){
    dots$fixed = start.value[names(dots$fixed)]
  }
  if (!"method" %in% names(dots)) dots$method = "Nelder-Mead"
  ctdf = aggregate(1:length(x),by=data.frame(abundance=x,gcn=gcn),FUN=length)
  #cat(sprintf("Fitting SAD on %d unique abundance/copy number pairs\n",dim(ctdf)[1]))
  if (is.null(trunc)) {
    LL = function(mu,sig,odisp){
      pps = dnblog(x = ctdf$abundance,gcn=ctdf$gcn, mu=mu, sig=exp(sig)+min.sig,
                   odisp=exp(odisp)+min.odisp, constant = constant,log = TRUE)*ctdf$x
      pp = -sum(pps)
      if(is.infinite(pp)) pp = .Machine$double.xmax
      return(pp)
    }
  }
  else {
    LL = function(mu,sig,odisp){
      pps = dtrunc_dnblog(x = ctdf$abundance,gcn=ctdf$gcn, mu=mu, sig=exp(sig)+min.sig,
                          odisp=exp(odisp)+min.odisp, trunc=trunc, constant = constant,log = TRUE)*ctdf$x
      pp = -sum(pps)
      if(is.infinite(pp)) pp = .Machine$double.xmax
      return(pp)
    }
  }
  result = do.call("mle2", c(list(LL, start = start.value,
                                  data = list(x = x)), dots))
  result@fullcoef["odisp"] = exp(result@fullcoef["odisp"])+min.odisp
  result@fullcoef["sig"] = exp(result@fullcoef["sig"])+min.sig
  return(result)
}
