#' @title  Poisson-lognormal distribution functions
#'
#' @description  Refined Poisson lognormal distribution functions for better accuracy after truncation
#'
#' @param x,q the observed species abundances; in 16S profiling, this is the observed 16S read count
#' @param n the number of random numbers to be generated
#' @param gcn the 16S copy number of OTUs; it should be a numeric vector of the same length as `x`
#' @param mu the log mean of the log-normal distribution
#' @param sig the log standard deviation of the log-normal distribution
#' @param trunc the truncation cutoff
#' @param constant the node and weight values for Gauss-Kronrod quadrature formula; see `Details` for more information
#' @param log,log.p logical; if `TRUE`, the log-transformed value will be returned; defaults to `FALSE`
#' @details `sads::dpoilog` encounters numerical instability that may lead to inaccurate likelihood for certain parameter combinations. The alternative functions implemented in this package aims to minimize the numerical error in model fitting.
#' @details `d_poilog` and `dtrunc_dpoilog` are the (truncated) probability mass function; `p_poilog` is the cumulative distribution function; `r_poilog` is the random number generator function.
#' @details `dconv_poilog` and  `pconv_poilog` are the integrand functions for `d_poilog` and `p_poilog`, respectively; despite their role as mostly internal functions, they are exported for user-created functions.
#' @details In Poisson-log-normal distribution, the sampling effort `eta`, has perfect collinearity with the `mu` in the log-normal distribution. Therefore, `eta` is degenerated and fixed to 1 internally (before GCN correction).
#' @details `constant` can be either a character value naming the specific element to use in the data list `quadgk_constant`, or a user-supplied list of numeric vectors containing the node and weight values.
#' @details The quantile function `q_poilog` is not implemented, and thus `r_poilog` can only generate random values using the full distribution without truncation, making it not suitable for RAD prediction.
#' @details For RAD prediction, use `representative_poilog` for approximated quantiles.
#' @details `fit_poilog` fits the Poisson-lognormal distribution using the refined PMF.
#' @return `d_poilog`, `dtrunc_dpoilog`, `p_poilog`, `dconv_poilog` and `pconv_poilog` return a numeric vector the same length as `x` or `q`
#' @return `r_poilog` return a numeric vector of length `n`
#' @return `fit_poilog` return a `mle2` class object of fitted model
#' @export
#' @rdname poilog
#'
d_poilog = function(x,gcn=rep(1,length(x)),mu,sig,constant="G25_K51",log=FALSE){
  pp = sapply(1:length(x),function(i){
    xx = x[i]
    if(xx<0) return(-Inf)
    if(is.infinite(xx)) return(-Inf)
    muc = mu+log(gcn[i])
    #
    peak.est = max(log(xx),muc)
    peak.det = function(y){y-(muc+(xx-1)*sig^2)+exp(y)*sig^2}
    try({peak.est = uniroot(f = peak.det,lower = min(muc,log(xx+1))-1,
                            upper = peak.est+1)$root},silent = TRUE)
    fin.intervals = peak.est+sig*seq(-5,5,1)
    fin.intervals = fin.intervals[is.finite(fin.intervals)]
    #
    ppp = integral_refined_02inf_lnorm(
      fun = dconv_poilog,params = list(q=xx,mu=muc,sig=sig),
      fin.intervals=fin.intervals,constant=constant,log=TRUE,silent = TRUE)
    return(ppp)
  })
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname poilog
#'
dconv_poilog = function(x,q,mu,sig,log=FALSE,adjust=0){
  pp =  dnorm(x = x,mean = mu,sd = sig,log = TRUE)+
    dpois(x = q,lambda = exp(x),log = TRUE)-
    adjust
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname poilog
#'
p_poilog = function(q,gcn=rep(1,length(q)),mu,sig,log.p=FALSE,lower.tail=FALSE,constant="G25_K51"){
  pp = sapply(1:length(q),function(i){
    qq = q[i]
    if(qq<0) return(0)
    if(is.infinite(qq)) return(-Inf)
    muc = mu+log(gcn[i])
    #
    peak.est = max(log(qq),muc)
    peak.det = function(y){y-(muc+qq*sig^2)+exp(y)*sig^2}
    try({peak.est = uniroot(f = peak.det,lower = min(muc,log(qq+1))-1,
                            upper = peak.est+1)$root},silent = TRUE)
    fin.intervals = peak.est+sig*seq(-5,5,1)
    fin.intervals = fin.intervals[is.finite(fin.intervals)]
    #
    ppp = integral_refined_02inf_lnorm(fun = pconv_poilog,params = list(q=qq,mu=muc,sig=sig),
                                       fin.intervals=fin.intervals,constant = constant,silent=TRUE,log = TRUE)
    return(ppp)
  })
  if(lower.tail) pp = log_one_minus_x(pp)
  if(log.p) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname poilog
#'
pconv_poilog = function(x,q,mu,sig,log=FALSE,adjust=0){
  pp =  dnorm(x = x,mean = mu,sd = sig,log = TRUE)+
    ppois(q = q,lambda = exp(x),log.p = TRUE,lower.tail = FALSE)-
    adjust
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname poilog
#'
dtrunc_dpoilog = function(x,gcn=rep(1,length(x)),mu,sig,trunc=0,constant="G25_K51",log=FALSE,debug=FALSE){
  if(min(x)<=trunc) stop("Truncation should be smaller than the minimum abundance!")
  unique.gcn = unique(gcn)
  gcn.idx = as.numeric(factor(gcn,levels = unique.gcn))
  unique.lower.sum = p_poilog(q = rep(trunc,length(unique.gcn)),gcn=unique.gcn,
                              mu = mu,sig = sig,
                              constant = constant,log.p = TRUE,lower.tail = FALSE)
  lower.sum = unique.lower.sum[gcn.idx]
  # initialize by the unobserved sum
  pp = -lower.sum
  # calculate truncated probability with exceptions when x==trunc+1
  if(any(x>(trunc+1))){
    these.x = x>(trunc+1)
    pp[these.x] = pp[these.x]+d_poilog(x = x[these.x],gcn=gcn[these.x],mu=mu,sig=sig,log=TRUE)
  }
  if(any(x==(trunc+1))){
    these.x = x==(trunc+1)
    pp[these.x] = log_one_minus_x(
      pp[these.x]+p_poilog(q = x[these.x],gcn=gcn[these.x],mu=mu,sig=sig,
                           constant=constant,lower.tail = FALSE,log.p=TRUE))
  }
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname poilog
#'
r_poilog = function(n,mu,sig){
  x = rlnorm(n = n,meanlog = mu,sdlog = sig)
  k = rpois(n = n,lambda = x)
  return(k)
}

#' @export
#' @rdname poilog
#'
fit_poilog = function(x,gcn=rep(1,length(x)),trunc=0,start.value,constant="G25_K51",min.sig = 1e-5,...){
  dots = list(...)
  if (any(x < 0) | any(!is.wholenumber(x)))
    stop("All x must be non-negative integers")
  if (!is.null(trunc)) {
    if (min(x) <= trunc)
      stop("truncation point should be lower than the lowest data value")
  }
  if (missing(start.value)){
    start.value = list(mu=0,sig=1)
  }
  start.value = start.value[c("mu","sig")]
  if ("fixed" %in% names(dots)){
    start.value[names(dots$fixed)] = dots$fixed
  }
  start.value$sig = log(start.value$sig)
  if ("fixed" %in% names(dots)){
    dots$fixed = start.value[names(dots$fixed)]
  }
  if (!"method" %in% names(dots)) dots$method = "Nelder-Mead"
  ctdf = aggregate(1:length(x),by=data.frame(abundance=x,gcn=gcn),FUN=length)
  #cat(sprintf("Fitting SAD on %d unique abundance/copy number pairs\n",dim(ctdf)[1]))
  if (is.null(trunc)) {
    LL = function(mu,sig){
      pps = d_poilog(x = ctdf$abundance,gcn=ctdf$gcn, mu=mu, sig=exp(sig)+min.sig,
                     constant = constant,log = TRUE)*ctdf$x
      pp = -sum(pps)
      if(is.infinite(pp)) pp = .Machine$double.xmax
      return(pp)
    }
  }
  else {
    LL = function(mu,sig){
      pps = dtrunc_dpoilog(x = ctdf$abundance,gcn=ctdf$gcn, mu=mu, sig=exp(sig)+min.sig,
                           trunc=trunc, constant = constant,log = TRUE)*ctdf$x
      pp = -sum(pps)
      if(is.infinite(pp)) pp = .Machine$double.xmax
      return(pp)
    }
  }
  result = do.call("mle2", c(list(LL, start = start.value,
                                  data = list(x = x)), dots))
  result@fullcoef["sig"] = exp(result@fullcoef["sig"])+min.sig
  return(result)
}
