#' @title  Probability functions of power-bend distribution
#'
#' @description  Alternative implementation of power-bend distribution functions with faster speed and better accuracy for large x or q
#'
#' @param x,q vector of quantiles (abundances or read numbers in a sample)
#' @param p vector of cumulative probabilities
#' @param n the number of random values to be generated
#' @param s exponent of the denominator in the power-bend distribution
#' @param lambda,minlogLambda the bending parameter of the distribution (scaling coefficient of the exponent in numerator)
#' @param adjust normalizing coefficient for integration of very large or small values
#' @param trunc truncation point of the distribution; defaults to `NULL`
#' @param exact the number of exact sum to be calculated in the infinite series
#' @param log,log.p logical;if the probabilities returned by or supplied to the function are in log-scale
#' @param lower.tail logical;if the cumulative probabilities are lower-tail (`TRUE`) or upper-tail (`FALSE`)
#' @param start.value the starting values of fitting a power-bend distribution
#' @param ... other parameters supplied to the `mle2` function
#' @details The `[dpqr]_powbend` functions are counterparts of the `[dpqr]powbend` functions in the `sads` package with better support and faster computation time for s<1 and large quantiles.
#' @details The `powbend_terms` function calculates the unnormalized terms in a power-bend distribution that are useful for many distribution functions based on power-bend distribution.
#' @details The `fit_powbend` function a counterpart of `sads::fitpowbend` that uses `d_powbend` function.
#' @return `powbend_term`,`d_powbend`,`p_powbend`,`r_powbend`,`q_powbend` returns a numeric vector
#' @return `fit_powbend` returns a object of `mle2` class
#' @export
#' @rdname d_powerbend
powbend_term = function(x,s,lambda,adjust=0,log=FALSE){
  dd = -lambda*x - s*log(x) + adjust
  if(log) return(dd)
  return(exp(dd))
}

#' @export
#' @rdname d_powerbend
d_powbend = function(x,s,lambda=0.001,minlogLambda=-log(lambda),trunc=NULL,exact=1e5,log=FALSE){
  if (!missing(lambda) && !missing(minlogLambda))
    stop("specify 'lambda' or 'minlogLambda' but not both")
  lambda = exp(-minlogLambda)
  if(is.null(trunc)){
    p0 = polylog_approximation(lambda = lambda,s = s,exact = exact,log = TRUE)
  }else{
    if(any(x<=trunc)) stop("x should be greater than trunc!")
    p0 = Lerch_transcendent_approximation(x=trunc+1,s=s,lambda=lambda,exact=exact,log=TRUE)
  }
  pp = powbend_term(x = x,s=s,lambda=lambda,log=TRUE)-p0
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname d_powerbend
p_powbend = function(q,s,lambda=0.001,minlogLambda=-log(lambda),trunc=NULL,exact=1e5,log.p=FALSE,lower.tail=TRUE){
  if (!missing(lambda) && !missing(minlogLambda))
    stop("specify 'lambda' or 'minlogLambda' but not both")
  lambda = exp(-minlogLambda)
  if(is.null(trunc)){
    p0 = polylog_approximation(s=s,lambda=lambda,exact=exact,log=TRUE)
    trunc=0
  }else{
    if(any(q<=trunc)) stop("q should be greater than trunc!")
    p0 = Lerch_transcendent_approximation(x=trunc+1,s=s,lambda=lambda,exact=exact,log=TRUE)
  }
  pp = sapply(q+1,function(qq){
    if(qq<1) return(p0)
    if(is.infinite(qq)) return(-Inf)
    Lerch_transcendent_approximation(x=qq,s=s,lambda=lambda,exact=exact,log=TRUE)
  })-p0
  pp.check = pp>0
  if(lower.tail){
    pp = log_one_minus_x(pp)
    if(any(pp.check)){
      pp[pp.check] = sapply(q[pp.check],function(qq){
        log_sum_x(powbend_term(x = (trunc+1):qq,s=s,lambda=lambda,log=TRUE))-p0
      })
    }
  }else{
    if(any(pp.check)){
      pp[pp.check] = sapply(q[pp.check],function(qq){
        log_one_minus_x(log_sum_x(powbend_term(x = (trunc+1):qq,s=s,lambda=lambda,log=TRUE))-p0)
      })
    }
  }
  if(log.p) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname d_powerbend
q_powbend = function(p,s,lambda=0.001,minlogLambda=-log(lambda),trunc=NULL,exact=1e5,log.p=FALSE,lower.tail=TRUE){
  if (!missing(lambda) && !missing(minlogLambda))
    stop("specify 'lambda' or 'minlogLambda' but not both")
  lambda = exp(-minlogLambda)
  z0 = polylog_approximation(lambda = lambda,s = s,exact = exact,log = TRUE)
  if(is.null(trunc)){
    p0 = z0
    trunc=0
  }else{
    p0 = Lerch_transcendent_approximation(x=trunc+1,s=s,lambda=lambda,exact=exact,log=TRUE)
  }
  if(!log.p) p = log(p)
  if(lower.tail) p = log_one_minus_x(p)
  p = sort(p,index.return=TRUE,decreasing = TRUE)
  p.idx = p$ix
  p = p$x + p0 - z0
  if(is.null(trunc)){
    qexact = 1:(exact+1)
  }else{
    qexact = 1+(trunc:exact)
  }
  if((is.wholenumber(s)&(s!=1))|(s<=0)){
    pexact = p_powbend(q = qexact,s = s,lambda = lambda,trunc = NULL,log.p = TRUE,lower.tail = FALSE)
  }else{
    pexact = ppowbend(q = qexact,s = s,omega = lambda,log.p = TRUE,lower.tail = FALSE)
  }
  if(any(is.na(pexact))) pexact[is.na(pexact)] = -Inf
  if(any(p>=min(pexact))){
    qq = sapply(p[p>=min(pexact)],function(pp){
      qexact[min(which(pp>=pexact))]
    })
  }else{
    qq = numeric(0)
  }
  if(length(qq)<length(p)){
    p2root = which(p<min(pexact))
    ff = function(xx,pp){
      pp - (Lerch_transcendent_approximation(x=xx,s=s,lambda=lambda,exact=exact,log=TRUE) - z0)
    }
    pexact = p_powbend(q = 10^seq(1,100,1),
                       s=s,lambda = lambda,trunc = NULL,exact = exact,log.p = TRUE,lower.tail = FALSE)
    qexact = c(10^seq(1,100,1),Inf)
    qq = c(qq,numeric(length(p)-length(qq)))
    for(i in p2root){
      if(all((p[i]-pexact)<0)){
        qq[i] = Inf
        next
      }
      qupper = qexact[min(which(p[i]>=pexact))]+1
      qlower = max(qq[i-1]-1,min(qexact))
      qq[i] = round(floor(uniroot(f = ff, pp = p[i], lower = qlower,upper = qupper,
                                  tol = sqrt(.Machine$double.eps),
                                  maxiter = 2000)$root))
    }
  }
  qq[p.idx] = qq
  return(qq)
}

#' @export
#' @rdname d_powerbend
r_powbend = function(n,s,lambda=0.001,minlogLambda=-log(lambda),exact=1e5,trunc=NULL){
  if (!missing(lambda) && !missing(minlogLambda))
    stop("specify 'lambda' or 'minlogLambda' but not both")
  lambda = exp(-minlogLambda)
  p0 = runif(n = n,min = 0,max = 1)
  x = q_powbend(p = p0,s = s,lambda = lambda,trunc = trunc,exact=exact,log.p = FALSE)
  return(x)
}

#' @export
#' @rdname d_powerbend
fit_powbend = function(x,trunc=NULL,exact=1e5,start.value,...){
  dots = list(...)
  if (any(x < 0) | any(!is.wholenumber(x)))
    stop("All x must be non-negative integers")
  if (!is.null(trunc)) {
    if (min(x) <= trunc)
      stop("truncation point should be lower than the lowest data value")
  }
  if (missing(start.value)){
    start.value = list(s=1,minlogLambda=6.9)
  }
  start.value = start.value[c("s","minlogLambda")]
  if ("fixed" %in% names(dots)){
    start.value[names(dots$fixed)] = dots$fixed
  }
  if ("fixed" %in% names(dots)){
    dots$fixed = start.value[names(dots$fixed)]
  }
  if (!"method" %in% names(dots)) dots$method = "Nelder-Mead"
  if (dots$method == "L-BFGS-B") {
    if (!"lower" %in% names(dots))
      dots$lower = c(s = log(-5), minlogLambda = 1)
    if (!"upper" %in% names(dots))
      dots$upper = c(s = log(5), minlogLambda = 16)
  }
  ctdf = aggregate(1:length(x),by=data.frame(abundance=x),FUN=length)
  #cat(sprintf("Fitting SAD on %d unique abundances\n",dim(ctdf)[1]))
  LL = function(s, minlogLambda){
    pps = d_powbend(x = ctdf$abundance,exact=exact,trunc = trunc,
                     s = s, minlogLambda = minlogLambda,log = TRUE)*ctdf$x
    pp = -sum(pps)
    if(is.infinite(pp)) pp = .Machine$double.xmax
    return(pp)
  }
  result = do.call("mle2", c(list(LL, start = start.value,
                                  data = list(x = x)), dots))
  return(result)
}
