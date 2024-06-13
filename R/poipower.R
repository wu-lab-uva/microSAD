#' @title  Power-law SAD with Poisson sampling error
#'
#' @description  The PMF, CDF, RNG and fitting function of power-law distributed communities with Poisson sampling error.
#'
#' @param x,q the observed species abundances, non-negative integers
#' @param p the cumulative probabilities of the species abundances, real values from 0 to 1
#' @param n the number of random values to be generated
#' @param gcn the 16S copy number of each species, with the same length and order of x, positive integers
#' @param eta the sampling effort, positive real number
#' @param trunc the truncation point to fit a distribution whose lower-tail is truncated; non-negative integer
#' @param s exponent of the denominator in power-law distribution; it is the same as parameter `s` in the `sads` package
#' @param exact the number of terms to be calculated exactly before approximation applies
#' @param constant the name of Gauss-Kronrod quadrature constants; see `Details`
#' @param log,log.p if log-transformed probability should be returned or supplied; default to `FALSE`
#' @param Nmax the maximum abundance beyond which the RAD is approximated by quantiles in a pareto distribution
#' @param lower.tail if the lower-tail cumulative probability should be returned or received; default to `FALSE` for accuracy of large quantiles
#' @param adjust log-adjustment value to avoid precision problem in integration
#' @param start.value a list of named numeric values as the starting values of ML fitting
#' @param upper the upper limit of parameters in the ML fitting
#' @param ... other parameters passed to the `mle2` function
#' @details `dpoipower` and `dtrunc_dpoipower` calculates the full and truncated PMF for Poisson-power-law distribution
#' @details `ppoipower` calculates the CDF of Poisson-powerlaw distribution. To achieve the best accuracy for SAD fitting, it caculates the upper tail of the CDF by default
#' @details `qpoipower` calculates the quantile function for the Poisson-power-law distribution. Given a probability value, it estimates the corresponding value of the distribution at that probability. Large quantiles are approximated using the continuous pareto distribution.
#' @details `rpoipower` generates random numbers under the Poisson-power-law distribution.
#' @details `constant` denotes the name of a set of values used by Gauss-Kronrod quadrature. This packages provides two sets of these values: `G25_K51`, the default; and `G7_K15`, the same as used in the GK method in `pracma::integral`
#' @details `pconv_poipower` is the integrand function used internally by `ppoipower`
#' @details `fit_poipower` fits the Poisson-power-law distribution to observed SAD using maximum likelihood. Correction for 16S rRNA gene copy number can be done by supplying the corresponding GCN to the `gcn` parameter.
#' @return `dpoipower`,`ppoipower` and `dtrunc_poipower` returns a numeric vector of the (log) probability
#' @return `qpoipower`,`rpoipower` returns a numeric vector of the (random) quantiles
#' @return `fitpoipower` returns a `mle2` class object containing the fitted model

#' @export
#' @rdname poipower
dtrunc_dpoipower = function(x, gcn=rep(1,length(x)), trunc=0, eta, s,
                            exact=1e5, constant="G25_K51", log=FALSE, debug=FALSE){
  if(eta < 0) stop("Poisson sampling effort should be non-negative!")

  if(min(x) <= trunc) stop("Truncation should be smaller than the minimum abundance!")
  x=round(x)
  # calculate probability mass not observed
  unique.gcn = unique(gcn)
  gcn.idx = as.numeric(factor(gcn, levels = unique.gcn))

  unique.lower.sum = ppoipower(q=rep(trunc, length(unique.gcn)),
                               eta=eta, s=s,
                               exact=exact, log=TRUE, lower.tail=FALSE)
  lower.sum = unique.lower.sum[gcn.idx]

  # calculate truncated probability with exceptions when x == trunc + 1
  pp = dpoipower(x = x, eta=eta, s=s,log=TRUE)-lower.sum
  if(any(x==(trunc+1))){
    these.x = which(x==(trunc+1))
    these.x = these.x[pp[these.x]>-log(2)]
    if(length(these.x)>0){
      pp[these.x] = log_one_minus_x(ppoipower(q = x[these.x],gcn=gcn[these.x],eta=eta,s=s,
                                              exact=exact,constant=constant,lower.tail = FALSE,log.p=TRUE)-
                                      lower.sum[these.x])
    }
  }
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname poipower
dpoipower = function(x,gcn=rep(1,length(x)), eta,s, exact=1e5, log=FALSE, debug=FALSE){
  if(s <= 1) return(rep(NaN,length(x)))
  x = round(x)
  pp = sapply(1:length(x), function(i){
    xx = x[i]
    if(is.infinite(xx)) return(-Inf)
    etac = eta*gcn[i]
    ppp = dpois(x = xx,lambda = etac,log = TRUE) +
      polylog_approximation(lambda=etac,s=s-xx,log=TRUE)
    return(ppp)
  }) - log(VGAM::zeta(s)) + eta
  if (log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname poipower
ppoipower = function(q, gcn=rep(1,length(q)), eta, s,
                     exact=1e5, constant="G25_K51",
                     log.p = FALSE, lower.tail=TRUE, debug=FALSE){
  if(s<=1) return(rep(NaN,length(q)))
  q = floor(q)
  zeta.coef = log(VGAM::zeta(s))
  pp = sapply(1:length(q), function(i){
    qq = q[i]
    if(qq<0) return(zeta.coef)
    if(is.infinite(qq)) return(-Inf)
    etac = eta*gcn[i]
    ppp.exact = log_sum_x(pconv_poipower(1:exact,q=qq,etac,s,log = TRUE))
    ppp.remainder = pconv_poipower(x = 1+exact,q=qq,s = s,eta = etac,log = TRUE)-log(2)
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
      fun=pconv_poipower, s = s, q = qq, eta = etac,
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
#' @rdname poipower
pconv_poipower=function(x,q, eta,s, log=FALSE, adjust = 0){
  pp =  ppois(q = q,lambda = x*eta,lower.tail=FALSE, log = TRUE) - s*log(x) - adjust
  if(log) return(pp)
  return(exp(pp))
}

#' @export
#' @rdname poipower
fit_poipower = function (x, gcn=rep(1,length(x)), trunc=0, start.value, upper = c(s=20), ...)
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
    etahat =0
  }
  else {
    shat = log(start.value$s -1 )
    etahat = log(start.value$eta)
  }
  ctdf = aggregate(1:length(x),by=data.frame(abundance=x, gcn= gcn),FUN=length)
  if (is.null(trunc)) {

    LL = function(s,eta) {
      s = 1+1e-10+exp(s)
      eta = exp(eta)
      pps = dpoipower(x=ctdf$abundance, gcn=ctdf$gcn, s=s, eta=eta, log = TRUE)*ctdf$x
      pp = -sum(pps)
      if(is.infinite(pp)) pp = .Machine$double.xmax
      return(pp)}
  } else{
    LL = function(s,eta) {

      s = 1+1e-10+exp(s)
      eta = exp(eta)
      pps = dtrunc_dpoipower(x=ctdf$abundance, gcn=ctdf$gcn, trunc = trunc, s=s, eta=eta, log = TRUE)*ctdf$x
      pp = -sum(pps)
      if(is.infinite(pp)) pp = .Machine$double.xmax
      return(pp)}
  }
  try({
    result = do.call("mle2", c(list(LL, start = list(s = shat, eta=etahat),
                                     data = list(x = x),method = 'L-BFGS-B' ),
                                dots))
    })
  if (!exists('result')){
    result = do.call("mle2", c(list(LL, start = list(s = shat, eta=etahat),
                                     data = list(x = x),method = 'Nelder-Mead'),
                                dots))
  } else if (result@details$convergence !=0) {
    cat("L-BFGS-B method didnt converge, Nelder-Mead tried instead\n")
    result = do.call("mle2", c(list(LL, start = list(s = shat, eta=etahat),
                                     data = list(x = x),method = 'Nelder-Mead'),
                                dots))
  }
  result@fullcoef["eta"] = exp(result@fullcoef["eta"])
  result@fullcoef["s"] = 1+1e-10+exp(result@fullcoef["s"])
  return(result)
}

#' @export
#' @rdname poipower
qpoipower = function(p, s, eta, trunc=0, log.p=FALSE,lower.tail=FALSE,
                     exact=1e5, constant="G25_K51", Nmax=.Machine$integer.max, debug = FALSE){
  # convert all desired CP to log-scale upper tail CP
  if(!log.p) p = log(p)
  if(lower.tail) p = log_one_minus_x(p)
  # sort desired CP from largest (q->0) to smallest (q->Inf)
  p.order = order(p,decreasing = TRUE,na.last = TRUE)
  # initialize output quantiles with input CP values to handle NA and NaN
  qq = p[p.order]
  # adjust desired CP by trunc
  if(!is.null(trunc)){
    p = p+ppoipower(q = trunc,eta = eta,s = s,
                    exact = exact,constant = constant,log.p = TRUE,lower.tail = FALSE)
  }
  # initialize guess for q
  if(!is.null(trunc)){
    guess.q = c(trunc,trunc+1,trunc+2)
  }else{
    guess.q = c(-1,0,1)
  }
  check.p = ppoipower(q = guess.q, eta = eta, s = s, exact = exact,
                      constant = constant,log.p = TRUE,lower.tail = FALSE)
  # start looping through CPs
  use.approx = FALSE
  for(i in p.order){
    desired.p = p[i]
    if(is.na(desired.p)){
      # terminate loop
      break
    }
    # check if the desired quantile can be found within the window
    check.val = c((check.p[2]<=desired.p)&(check.p[1]>desired.p),
                  (check.p[3]<=desired.p)&(check.p[2]>desired.p),
                  check.p[3]>desired.p)
    if(check.val[3]){
      # the guessed window should be moved by >1
      # initialize the sampling points
      if(guess.q[3]<Nmax){
        minmax.q = c(guess.q[3],
                     min(Nmax,max(guess.q[3]+1,
                                  ceiling(qpareto(p = desired.p-check.p[3]+
                                                    ppareto(q=guess.q[3], shape = s-1,scale = 1,
                                                            lower.tail = FALSE,log.p = TRUE),
                                                  shape = s-1,scale = 1,
                                                  lower.tail = FALSE,log.p = TRUE)))))
      }else{
        # terminate loop and use vectorized approximation by qpareto
        use.approx = TRUE
        break
      }
      # calculates upper tail CP for them
      minmax.p = c(check.p[3],
                   ppoipower(q = minmax.q[2],eta = eta,s = s,exact = exact,
                             constant = constant,log.p = TRUE,lower.tail = FALSE))
      # make a new guess based on the sampled upper tail CP
      # ceiling function guarantees that the new guess is at least +2
      if((minmax.p[2]-minmax.p[1])==0){
        guess.q = guess.q+2
      }else{
        guess.q = ceiling(minmax.q[1] + (minmax.q[2]-minmax.q[1])*
                            ((desired.p-minmax.p[1])/(minmax.p[2]-minmax.p[1])))+c(-1,0,1)
      }
      if(guess.q[3]>Nmax){
        # limit the new guess to Nmax
        guess.q = ceiling(Nmax)+c(-1,0,1)
      }
      check.p = ppoipower(q = guess.q,eta = eta,s = s,exact = exact,
                          constant = constant,log.p = TRUE,lower.tail = FALSE)
      check.val = c(check.p[1]<=desired.p,
                    check.p[3]>desired.p)
      if((guess.q[3]>=Nmax)&(check.val[2])){
        # terminate loop and use vectorized approximation by qpareto
        use.approx=TRUE
        break
      }
      while(any(check.val)){
        # if the new guess is still off, update the sample points
        # and updates upper tail CP for them
        minmax.q = minmax.q[which.min(abs(minmax.p-desired.p))]
        minmax.p = minmax.p[which.min(abs(minmax.p-desired.p))]
        if(minmax.q==guess.q[2]){
          minmax.q = guess.q[-1]
          minmax.p = check.p[-1]
        }else{
          minmax.q = c(minmax.q,guess.q[2])
          minmax.p = c(minmax.p,
                       check.p[2])
        }
        # make a new guess again
        if((minmax.p[2]-minmax.p[1])==0){
          guess.q = guess.q+sum(as.numeric(check.val)*c(-1,1))
        }else{
          guess.q = ceiling(minmax.q[1] + (minmax.q[2]-minmax.q[1])*
                              ((desired.p-minmax.p[1])/(minmax.p[2]-minmax.p[1])))+c(-1,0,1)
        }
        # update the guessed window
        check.p = ppoipower(q = guess.q, eta = eta, s = s, exact = exact,
                            constant = constant,log.p = TRUE,lower.tail = FALSE)
        check.val = c(check.p[1]<=desired.p,
                      check.p[3]>desired.p)
      }
      # after a new window is reached, check the interval where the desired quantile lies
      check.val = c((check.p[2]<=desired.p)&(check.p[1]>desired.p),
                    (check.p[3]<=desired.p)&(check.p[2]>desired.p))
    }
    if(check.val[2]){
      # the guessed window moves by 1
      qq[i] = guess.q[3]
      guess.q = guess.q+1
      check.p = c(check.p[-1],
                  ppoipower(q = guess.q[3],eta = eta,s = s,exact = exact,
                            constant = constant,log.p = TRUE,lower.tail = FALSE))
    }else{
      # the guessed window stays
      qq[i] = guess.q[2]
    }
  }
  if(use.approx){
    # approximation by qpareto
    qq[p.order[-(1:(which(p.order==i)-1))]] =
      ceiling(qpareto(p = p[p.order[-(1:(which(p.order==i)-1))]]-check.p[3]+
                        ppareto(q=guess.q[3], shape = s-1,scale = 1,lower.tail = FALSE,log.p = TRUE),
                      shape = s-1,scale = 1,lower.tail = FALSE,log.p = TRUE))
  }
  return(qq)
}

#' @export
#' @rdname poipower
rpoipower = function(n, s, eta, trunc=0,
                     exact=1e5, constant="G25_K51", Nmax=.Machine$integer.max, debug = FALSE){
  pp = runif(n = n,min = 0,max = 1)
  qq = qpoipower(p = pp,s = s,eta = eta,trunc = trunc,log.p = FALSE,lower.tail = FALSE,
                 exact = exact,constant = constant,Nmax = Nmax,debug = debug)
  return(qq)
}