#' @title  Convert SAD fit to octave fit
#'
#' @description  Convert the SAD fit to octave fit to allow comparison with octave-based models like gambin
#'
#' @param x the species abundance for each species 
#' @param model the base name of model (e.g., poipowbend for Poisson-powerbend)
#' @param params the fitted full parameters of the model (e.g., from the `fullcoef` slot of a `mle2` class object)
#' @param pred_rad the predicted species abundance using the fitted model (e.g., from representative_* functions)
#' @param fixed the number of fixed parameters for calculating AIC and BIC
#' @details The octaves start with '0', which contains all the species with 1 individuals. Octave 'i' contains species with [2^i,2^(i+1)-1] individuals.
#' @details For SAD in octaves, the truncation is fixed to be 0.
#' @return `sad2octave` returns a list containing the SAD in octave representation and the corresponding fitting metrics if a SAD fit is provided
#' @export
#' @rdname sad2octave
sad2octave = function(x,model,params,pred_rad,fixed=0){
  # convert SAD to octaves
  octaves = seq(0,floor(log2(max(x))),1)
  octave.x = sapply(octaves,function(i){
    sum((x>=2^i)&(x<2^(i+1)))
  })
  names(octave.x) = octaves
  # return if no SAD fit is supplied
  if(missing(model)|missing(params)) return(list(octave=octave.x))
  # convert SAD fit to octave fit
  pfunc = paste0("p",model)
  dfunc = paste0("d",model)
  p0 = do.call(what = pfunc,args = c(list(q=0,lower.tail=FALSE,log.p=TRUE),params))
  octave.p = sapply(octaves,function(i){
    if(i<5){
      minx = 2^i
      maxx = 2^(i+1)-1
      pp = log_sum_x(do.call(what = dfunc,args = c(list(x=minx:maxx,log=TRUE),params)))
    }else{
      minx = 2^i
      maxx = 2^(i+1)-1
      pp = do.call(what = pfunc,args = c(list(q=c(minx-1,maxx-1),lower.tail=FALSE,log.p=TRUE),params))
      return(pp[1]+log_one_minus_x(pp[2]-pp[1]))
    }
  })-p0
  minuslogl = -sum(octave.p*octave.x)
  octave.AIC = 2*minuslogl+2*(length(params)-fixed)
  octave.BIC = 2*minuslogl+log(length(x))*(length(params)-fixed)
  octave.r2m = r2modified(obs = octave.x/sum(octave.x),pred = exp(octave.p),log = FALSE)
  if(missing(pred_rad)){
    octave.r2m_rad = NA
  }else{
    octave.rad = do.call(c,mapply(rep,octaves,octave.x,SIMPLIFY = FALSE))
    octave.rad_pred = floor(log2(sort(pred_rad)))
    octave.r2m_rad = r2modified(obs = octave.rad,pred = octave.rad_pred,log = FALSE)
  }
  return(list(octave=octave.x,probs=octave.p,minuslogl=minuslogl,
              AIC=octave.AIC,BIC=octave.BIC,r2m=octave.r2m,r2m_rad=octave.r2m_rad))
}

#' @title  Goodness-of-fit for octave fit by RAD-RAD comparison 
#'
#' @description  Check goodness-of-fit for octave fit in RAD-RAD comparison format to allow comparison with sad-based models
#'
#' @param octave.x the number of species in each octave (e.g., $Data$species from the gambin fit object)
#' @param octaves the octave index matching the `octave.x` argument (e.g., $Data$octave from the gambin fit object)
#' @param fitted.values the expected proportion of species in each octave, calculated from the model fit (e.g., $fitted.values from the gambin fit object)
#' @details The octaves start with '0', which contains all the species with 1 individuals. Octave 'i' contains species with [2^i,2^(i+1)-1] individuals.
#' @details `fitted.values` generate the empirical cumulative distribution function from which an RAD can be sampled like in representative_* functions.
#' @return `r2m_rad` returns a numeric value
#' @export
#' @rdname r2m_rad
r2m_rad = function(octave.x,octaves,fitted.values) {
  obs.octave.rad=do.call(
    c,
    mapply(FUN = function(x,y){rep(y,x)},
           octave.x,octaves,SIMPLIFY = FALSE))
  fitted.values = fitted.values/sum(fitted.values)
  this.ecdf = c(0,cumsum(fitted.values),Inf)
  this.ecdf.octave = c(0,octaves,max(octaves))
  rep.pp =(1:sum(octave.x)-0.5)/sum(octave.x)
  rep.octave.rad = sapply(rep.pp,function(pp){
    this.ecdf.octave[min(which(this.ecdf>=pp))]
  })
  this.r2m = r2modified(obs = sort(obs.octave.rad),
                        pred = sort(rep.octave.rad),log = FALSE)
  return(this.r2m)
}
