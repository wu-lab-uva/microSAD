#' @title  Convert SAD fit to octave fit
#'
#' @description  Convert the SAD fit to octave fit to allow comparison with octave-based models like gambin
#'
#' @param x the species abundance for each species 
#' @param model the predicted values
#' @param params the weight of each value; defaults to equal weights; these weights are used to balance data points from datasets of different sizes.
#' @param fixed the number of fixed parameters for calculating AIC and BIC
#' @details The octaves start with '0', which contains all the species with 1 individuals. Octave 'i' contains species with [2^i,2^(i+1)-1] individuals.
#' @details For SAD in octaves, the truncation is fixed to be 0.
#' @return `sad2octave` returns a list containing the SAD in octave representation and the corresponding fitting metrics if a SAD fit is provided
#' @export
#' @rdname r2modified
sad2octave = function(x,model,params,fixed=0){
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
  return(list(octave=octave.x,probs=octave.p,minuslogl=minuslogl,
              AIC=octave.AIC,BIC=octave.BIC,r2m=octave.r2m))
}

