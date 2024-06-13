#' @title  Generate RAD from SAD model
#'
#' @description  A set of functions for fast RAD approximation from SAD models
#'
#' @param n the number of observed species to sample from SAD
#' @param trunc the truncation point only above which the species are observed
#' @param s the exponent of the denominator x^s
#' @param lambda,minlogLambda the coefficient of the exponent of the numerator e^(-lambda*x)
#' @param alpha,N the alternative parameterization of log-series that is compatible with `sads::fitls`
#' @param oM alternative parameterization of `minlogLambda` that is compatible with `sads::fitpowbend`
#' @param eta the coefficient of sampling effort as the mean of the compounding Poisson/negative-binomial distribution
#' @param odisp the coefficient of sampling bias as the overdispersion in negative-binomial distribution
#' @param random logical; when `TRUE`, the first RAD is randomly (instead of evenly) sampled from the SAD
#' @param Nrep the number of replicates to be generated
#' @param Nexact the maximum abundance beyond which the RAD is approximated
#' @param exact the number of terms to be summed in the infinite series before approximation applies
#' @param maxit the maximum iteration before bounded quantile search is terminated
#' @param constant a character name for the Gauss-Kronrod system or a list of user-supplied numeric vectors; see `integral_refined()` for more information
#' @details `representative_*` are a set of functions dedicated to generating RAD replicates from SAD models to evaluate the absolute goodness of fit
#' @details They apply a hybrid strategy of exact and approximate quantile to balance the accuracy and time complexity of prediction. Specifically, larger quantiles are approximated by a simpler distribution for faster prediction without losing too much accuracy.
#' @return `representative_*` functions return a list of RAD whose elements are sorted numeric values
#' @export
#' @rdname representative_RAD
representative_power = function(n,trunc=NULL,s,random=FALSE,
                                Nrep =1, Nexact=1e6){
  all.rad = lapply(1:Nrep,function(ii){
    # desired cumulative probability for n quantiles
    if(random|ii>1){
      pp0 = sort(runif(n = n,min = 0,max = 1),decreasing = FALSE)
    }else{
      pp0 = ((1:n)-0.5)/n
    }
    pp0 = log(pp0)
    # adjust cumulative probability by truncation
    if(!is.null(trunc)){
      ppt = ppower(q=trunc,s=s,lower.tail=FALSE,log.p=TRUE)
      pp0 = pp0+ppt
    }
    # maximum cumulative probability and quantile before approximation applies
    pexact = ppower(q = Nexact,s = s,lower.tail = FALSE,log.p = TRUE)
    # calculate exact quantiles
    if(any(pp0>=pexact)){
      xexact = qpower(p = exp(pp0[pp0>=pexact]),s = s,lower.tail = FALSE,log.p = FALSE)
    }else{
      xexact = numeric(0)
    }
    # check the number of approximated quantiles
    napprox = n-length(xexact)
    if(napprox<=0) return(xexact)
    # adjust cumulative probability for approximation
    pp1 = pp0[pp0<pexact]
    pp1 = pp1-pexact+ppareto(q = Nexact,shape = s-1,scale = 1,lower.tail = FALSE,log.p = TRUE)
    # use the continuous Pareto distribution to approximate power distribution
    xapprox = round(qpareto(p = pp1, scale = 1,shape = s-1,lower.tail = FALSE,log.p = TRUE))
    return(sort(c(xexact,xapprox)))
  })
  return(all.rad)
}

#' @export
#' @rdname representative_RAD
representative_poipower = function(n,trunc=0,s,eta,random=FALSE,
                                Nrep =1, Nexact=1e4,constant="G25_K51"){
  # build exact reference CDF that takes time to compute
  if(!is.null(trunc)){
    qexact = trunc:Nexact
    pexact = ppoipower(q = qexact, s=s, eta = eta, log.p = TRUE,
                    lower.tail = FALSE, constant=constant)
    ppt = pexact[1]
    pexact = pexact-ppt
    qexact[1] = trunc+1
  }else{
    qexact = 0:Nexact
    pexact = ppoipower(q = qexact, s=s, eta = eta, log.p = TRUE,
                       lower.tail = FALSE, constant=constant)
    ppt = 0
    pexact = c(0,pexact)
    qexact = c(0,qexact)
  }
  #
  all.rad = lapply(1:Nrep,function(ii){
    # desired cumulative probability for n quantiles
    if(random|ii>1){
      pp0 = sort(runif(n = n,min = 0,max = 1),decreasing = FALSE)
    }else{
      pp0 = ((1:n)-0.5)/n
    }
    pp0 = log(pp0)
    if(any(pp0>=min(pexact))){
      xexact = sapply(pp0[pp0>=min(pexact)],function(ppp){
        qexact[min(which(ppp>=pexact))]
      })
    }else{
      xexact = numeric(0)
    }
    napprox = n-length(xexact)
    if(napprox<=0) return(rev(xexact))
    #
    pp1 = pp0[pp0<min(pexact)]
    pp1 = pp1-min(pexact)+
      ppareto(q = Nexact,shape = s-1,scale = 1,lower.tail = FALSE,log.p = TRUE)
    xapprox = round(qpareto(p = pp1,shape = s-1,scale = 1,lower.tail = FALSE,log.p = TRUE))
    return(sort(rev(c(xapprox,xexact))))
  })
  return(all.rad)
}

#' @export
#' @rdname representative_RAD
representative_ls = function(n,trunc=NULL,lambda=0.01,minlogLambda=-log(lambda),
                             alpha, N,
                             s, oM,
                             random=FALSE, Nrep = 1, Nexact=1e4){
  if (!missing(lambda) && !missing(minlogLambda))
    stop("specify 'lambda' or 'minlogLambda' but not both")
  lambda = exp(-minlogLambda)
  # handling output from sads::fitls
  if(!missing(alpha) && !missing(N)){
    lambda = -log(N/(alpha+N))
  }
  # handling output from sads::fitpowbend
  if(!missing(s) && !missing(oM)){
    lambda = exp(-oM)
    if(s!=1) warning("Detecting s!=1 in log-series, forcing s to 1.",immediate. = TRUE)
  }
  return(representative_powbend(n=n,trunc=trunc,lambda=lambda,s=1,
                                random=random,Nrep=Nrep,Nexact=Nexact))
}

#' @export
#' @rdname representative_RAD
representative_powbend = function(n,trunc=NULL,s,lambda=0.01,minlogLambda=-log(lambda),
                                  omega, oM,
                                  random=FALSE, Nrep = 1, Nexact=1e4){
  if (!missing(lambda) && !missing(minlogLambda))
    stop("specify 'lambda' or 'minlogLambda' but not both")
  if (!missing(omega) && !missing(oM))
    stop("specify 'omega' or 'oM' but not both")
  if(!missing(oM)) omega = exp(-oM)
  if (!missing(omega) || !missing(oM)){
    lambda = omega
  }else{
    lambda = exp(-minlogLambda)
  }
  all.rad = lapply(1:Nrep,function(ii){
    if(random|ii>1){
      pp0 = sort(runif(n = n,min = 0,max = 1),decreasing = FALSE)
    }else{
      pp0 = ((1:n)-0.5)/n
    }
    xx = q_powbend(p = pp0,s = s,lambda = lambda,trunc = trunc,exact = Nexact,log.p = FALSE,lower.tail = FALSE)
    return(rev(xx))
  })
  return(all.rad)
}

#' @export
#' @rdname representative_RAD
representative_poilog = function(n,trunc=0,mu=0,sig=1,random=FALSE,
                                 Nrep = 1, Nexact=1e4,constant="G25_K51"){
  # build exact reference probability that takes time to compute
  if(!is.null(trunc)){
    qexact = trunc:Nexact
    pexact = p_poilog(q = qexact, mu = mu, sig = sig, log.p = TRUE,
                      lower.tail = FALSE, constant = constant)
    ppt = pexact[1]
    pexact = pexact-ppt
    qexact[1] = trunc+1
  }else{
    qexact = 0:Nexact
    pexact = p_poilog(q = qexact, mu = mu, sig = sig, log.p = TRUE,
                      lower.tail = FALSE, constant = constant)
    ppt = 0
    pexact = c(0,pexact)
    qexact = c(0,qexact)
  }
  #
  all.rad = lapply(1:Nrep,function(ii){
    if(random|ii>1){
      pp0 = sort(runif(n = n,min = 0,max = 1),decreasing = FALSE)
    }else{
      pp0 = ((1:n)-0.5)/n
    }
    pp0 = log(pp0)
    if(any(pp0>=min(pexact))){
      xexact = sapply(pp0[pp0>=min(pexact)],function(ppp){
        qexact[min(which(ppp>=pexact))]
      })
    }else{
      xexact = numeric(0)
    }
    napprox = n-length(xexact)
    if(napprox<=0) return(rev(xexact))
    #
    pp1 = pp0[pp0<min(pexact)]
    pp1 = pp1-min(pexact)+
      plnorm(q = Nexact,meanlog = mu,sdlog = sig,lower.tail = FALSE,log.p = TRUE)
    xapprox = round(qlnorm(p = pp1,meanlog = mu,sdlog = sig,lower.tail = FALSE,log.p = TRUE))
    return(sort(rev(c(xapprox,xexact))))
  })
  return(all.rad)
}

#' @export
#' @rdname representative_RAD
representative_nblog = function(n,trunc=0,mu=0,sig=1,odisp=1,random=FALSE,
                                Nrep = 1, Nexact=1e4, maxit=16, constant="G25_K51"){
  # build exact reference CDF that takes time to compute
  if(!is.null(trunc)){
    qexact = trunc:Nexact
    pexact = pnblog(q = qexact, mu = mu, sig = sig, odisp = odisp, log.p = TRUE,
                    lower.tail = FALSE, constant=constant)
    ppt = pexact[1]
    pexact = pexact-ppt
    qexact[1] = trunc+1
  }else{
    qexact = 0:Nexact
    pexact = pnblog(q = qexact, mu = mu, sig = sig, odisp = odisp, log.p = TRUE,
                    lower.tail = FALSE, constant=constant)
    ppt = 0
    pexact = c(0,pexact)
    qexact = c(0,qexact)
  }
  #
  all.rad = lapply(1:Nrep,function(ii){
    if(random|ii>1){
      pp0 = sort(runif(n = n,min = 0,max = 1),decreasing = FALSE)
    }else{
      pp0 = ((1:n)-0.5)/n
    }
    pp0 = log(pp0)
    if(any(pp0>=min(pexact))){
      xexact = sapply(pp0[pp0>=min(pexact)],function(ppp){
        qexact[min(which(ppp>=pexact))]
      })
    }else{
      xexact = numeric(0)
    }
    napprox = n-length(xexact)
    if(napprox<=0) return(rev(xexact))
    #
    pp1 = pp0[pp0<min(pexact)]
    pp2 = pp1-min(pexact)+
      plnorm(q = Nexact+1,meanlog = mu,sdlog = sig,lower.tail = FALSE,log.p = TRUE)
    xapprox1 = round(qlnorm(p = pp2,meanlog = mu,sdlog = sig,lower.tail = FALSE,log.p = TRUE))
    papprox1 = pnblog(q = c(xapprox1,2*max(xapprox1)), mu = mu, sig = sig, odisp = odisp, log.p = TRUE,
                      lower.tail = FALSE, constant=constant)-ppt
    qexact1 = c(max(qexact),rev(xapprox1),2*max(xapprox1))
    pexact1 = c(min(pexact),rev(papprox1)[-1],min(papprox1))
    while((min(pexact1)>min(pp1))){
      qexact1 = c(qexact1,2*max(qexact1))
      pexact1 = c(pexact1,
                  pnblog(q = 2*max(qexact1), mu = mu, sig = sig, odisp = odisp, log.p = TRUE,
                         lower.tail = FALSE, constant=constant)-ppt)
    }
    try.count = 0
    papprox1 = papprox1[1:length(pp1)]
    while((max(abs(pp1-papprox1))>1e-6)&(try.count<=maxit)){
      for(i in 1:length(pp1)){
        leftapprox.idx = max(which(pexact1>pp1[i]))
        rightapprox.idx = min(which(pexact1<=pp1[i]))
        xapprox1[i] = round((qexact1[leftapprox.idx]*(pp1[i]-pexact1[rightapprox.idx])
                             +qexact1[rightapprox.idx]*(pexact1[leftapprox.idx]-pp1[i]))/
                              (pexact1[leftapprox.idx]-pexact1[rightapprox.idx]))
        if(xapprox1[i]%in%qexact1){
          papprox1[i] = pexact1[which(qexact1==xapprox1[i])]
        }else{
          papprox1[i] = pnblog(q = xapprox1[i], mu = mu, sig = sig, odisp = odisp, log.p = TRUE,
                               lower.tail = FALSE, constant=constant)-ppt
          qexact1 = sort(c(qexact1,xapprox1[i]),decreasing = FALSE)
          pexact1 = sort(c(pexact1,papprox1[i]),decreasing = TRUE)
        }
      }
      try.count = try.count+1
    }
    return(sort(rev(c(xapprox1,xexact))))
  })
  return(all.rad)
}

#' @export
#' @rdname representative_RAD
representative_poipowbend = function(n,trunc=0,s,lambda=0.01,minlogLambda=-log(lambda),eta=1e-3,
                                     random=FALSE, Nrep = 1, Nexact = 1e4,
                                     exact=1e5, constant = "G25_K51"){
  if(eta<0) stop("Poisson sampling effort should be non-negative!")
  if (!missing(lambda) && !missing(minlogLambda))
    stop("specify 'lambda' or 'minlogLambda' but not both")
  lambda = exp(-minlogLambda)
  # build exact reference probability that takes time to compute
  if(!is.null(trunc)){
    qexact = trunc:Nexact
    pexact = ppoipowbend(q = qexact, s = s, lambda = lambda, eta = eta,
                         exact = exact, constant=constant, log.p = TRUE, lower.tail = FALSE)
    ppt = pexact[1]
    pexact = pexact-ppt
    qexact[1] = trunc+1
  }else{
    qexact = 0:Nexact
    pexact = ppoipowbend(q = qexact, s = s, lambda = lambda, eta = eta,
                         exact = exact, constant=constant, log.p = TRUE, lower.tail = FALSE)
    ppt = 0
    pexact = c(0,pexact)
    qexact = c(0,qexact)
  }
  #
  all.rad = lapply(1:Nrep,function(ii){
    if(random|ii>1){
      pp0 = sort(runif(n = n,min = 0,max = 1),decreasing = FALSE)
    }else{
      pp0 = ((1:n)-0.5)/n
    }
    pp0 = log(pp0)
    if(any(pp0>=min(pexact))){
      xexact = sapply(pp0[pp0>=min(pexact)],function(ppp){
        qexact[min(which(ppp>=pexact))]
      })
    }else{
      xexact = numeric(0)
    }
    napprox = n-length(xexact)
    if(napprox<=0) return(rev(xexact))
    #
    pp1 = pp0[pp0<min(pexact)]
    if((lambda+eta)<0.9){
      pp1 = pp1-min(pexact)+
        p_powbend(q = Nexact,s = s,lambda = log1p(lambda/eta),
                  exact = exact,log.p = TRUE,lower.tail = FALSE)
      xapprox = q_powbend(p = pp1, s = s, lambda = log1p(lambda/eta), log.p = TRUE,lower.tail = FALSE)
    }else{
      pp1 = pp1-min(pexact)+
        p_powbend(q = Nexact,s = s,lambda = lambda/eta,
                  exact = exact,log.p = TRUE,lower.tail = FALSE)
      xapprox = q_powbend(p = pp1, s = s, lambda = lambda/eta, log.p = TRUE,lower.tail = FALSE)
    }
    return(sort(rev(c(xapprox,xexact))))
  })
  return(all.rad)
}

#' @export
#' @rdname representative_RAD
representative_nbpowbend = function(n,trunc=0,s,lambda=0.01,minlogLambda=-log(lambda),eta=1e-3,odisp=1,
                                    random=FALSE, Nrep = 1, Nexact = 1e4, maxit=16,
                                    exact=1e5, constant = "G25_K51"){
  if(eta<0) stop("Poisson sampling effort should be non-negative!")
  if (!missing(lambda) && !missing(minlogLambda))
    stop("specify 'lambda' or 'minlogLambda' but not both")
  lambda = exp(-minlogLambda)
  # build exact reference CDF that takes time to compute
  if(!is.null(trunc)){
    qexact = trunc:Nexact
    pexact = pnbpowbend(q = qexact, s = s, lambda = lambda, eta = eta, odisp = odisp,
                         exact = exact, constant=constant, log.p = TRUE, lower.tail = FALSE)
    ppt = pexact[1]
    pexact = pexact-ppt
    qexact[1] = trunc+1
  }else{
    qexact = 0:Nexact
    pexact = pnbpowbend(q = qexact, s = s, lambda = lambda, eta = eta, odisp = odisp,
                         exact = exact, constant=constant, log.p = TRUE, lower.tail = FALSE)
    ppt = 0
    pexact = c(0,pexact)
    qexact = c(0,qexact)
  }
  #
  all.rad = lapply(1:Nrep,function(ii){
    if(random|ii>1){
      pp0 = sort(runif(n = n,min = 0,max = 1),decreasing = FALSE)
    }else{
      pp0 = ((1:n)-0.5)/n
    }
    pp0 = log(pp0)
    if(any(pp0>=min(pexact))){
      xexact = sapply(pp0[pp0>=min(pexact)],function(ppp){
        qexact[min(which(ppp>=pexact))]
      })
    }else{
      xexact = numeric(0)
    }
    napprox = n-length(xexact)
    if(napprox<=0) return(rev(xexact))
    #
    pp1 = pp0[pp0<min(pexact)]
    if((lambda+eta)<0.9){
      pp2 = pp1-min(pexact)+
        p_powbend(q = Nexact,s = s,lambda = log1p(lambda/eta),
                  exact = exact,log.p = TRUE,lower.tail = FALSE)
      xapprox1 = q_powbend(p = pp2, s = s, lambda = log1p(lambda/eta), log.p = TRUE,lower.tail = FALSE)
    }else{
      pp2 = pp1-min(pexact)+
        p_powbend(q = Nexact,s = s,lambda = lambda/eta,
                  exact = exact,log.p = TRUE,lower.tail = FALSE)
      xapprox1 = q_powbend(p = pp2, s = s, lambda = lambda/eta, log.p = TRUE,lower.tail = FALSE)
    }
    papprox1 = pnbpowbend(q = c(xapprox1,2*max(xapprox1)), s = s, lambda = lambda, eta = eta, odisp = odisp,
                          exact = exact, constant=constant, log.p = TRUE, lower.tail = FALSE)-ppt
    qexact1 = c(max(qexact),rev(xapprox1),2*max(xapprox1))
    pexact1 = c(min(pexact),rev(papprox1)[-1],min(papprox1))
    while((min(pexact1)>min(pp1))){
      qexact1 = c(qexact1,2*max(qexact1))
      pexact1 = c(pexact1,
                  pnbpowbend(q = 2*max(qexact1), s = s, lambda = lambda, eta = eta, odisp = odisp,
                             exact = exact, constant=constant, log.p = TRUE, lower.tail = FALSE)-ppt)
    }
    try.count = 0
    papprox1 = papprox1[1:length(pp1)]
    while((max(abs(pp1-papprox1))>1e-6)&(try.count<=maxit)){
      for(i in 1:length(pp1)){
        leftapprox.idx = max(which(pexact1>pp1[i]))
        rightapprox.idx = min(which(pexact1<=pp1[i]))
        xapprox1[i] = round((qexact1[leftapprox.idx]*(pp1[i]-pexact1[rightapprox.idx])
                             +qexact1[rightapprox.idx]*(pexact1[leftapprox.idx]-pp1[i]))/
                              (pexact1[leftapprox.idx]-pexact1[rightapprox.idx]))
        if(xapprox1[i]%in%qexact1){
          papprox1[i] = pexact1[which(qexact1==xapprox1[i])]
        }else{
          papprox1[i] = pnbpowbend(q = xapprox1[i], s = s, lambda = lambda, eta = eta, odisp = odisp,
                                   exact = exact, constant=constant, log.p = TRUE, lower.tail = FALSE)-ppt
          qexact1 = sort(c(qexact1,xapprox1[i]),decreasing = FALSE)
          pexact1 = sort(c(pexact1,papprox1[i]),decreasing = TRUE)
        }
      }
      try.count = try.count+1
    }
    return(sort(rev(c(xapprox1,xexact))))
  })
  return(all.rad)
}

#' @export
#' @rdname representative_RAD
representative_poils = function(...){
  dots = list(...)
  if(is.null(dots$s)) dots$s = 1
  res = do.call("representative_poipowbend",dots)
  return(res)
}

#' @export
#' @rdname representative_RAD
representative_nbls = function(...){
  dots = list(...)
  if(is.null(dots$s)) dots$s = 1
  res = do.call("representative_nbpowbend",dots)
  return(res)
}

