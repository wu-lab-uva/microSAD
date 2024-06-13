#' @title  Goodness of fit by R^2
#'
#' @description  Calculate the R^2 of a 1 to 1 fit
#'
#' @param obs the observed values
#' @param pred the predicted values
#' @param weight the weight of each value; defaults to equal weights; these weights are used to balance data points from datasets of different sizes.
#' @param log if `TRUE`, the R^2 is calculated on a log-log scale; defaults to `FALSE`.
#' @param normalize if `TRUE`, the R^2 is normalized between 0 and 1; defaults to `FALSE`.
#' @details This function is adapted from the modified coefficient of determination described by Shoemaker et al (2017). Instead of coefficient of determination R^2 in linear regression, `r2modified` calculates the R^2 of 1 to 1 fit; as a result, it can be negative.
#' @references Shoemaker, William R., Kenneth J. Locey, and Jay T. Lennon. "A macroecological theory of microbial biodiversity." Nature ecology & evolution 1.5 (2017): 0107.
#' @return `r2modified` returns a single numeric value <=1
#' @export
#' @rdname r2modified
r2modified = function(obs,pred,weight=rep(1,length(obs)),log=FALSE,normalize=FALSE){
  if(log){
    obs = log10(obs)
    pred = log10(pred)
  }
  mean.obs = sum(obs*weight)/sum(weight)
  r2m = sum(((obs-pred)^2)*weight)/sum(((obs-mean.obs)^2)*weight)
  if(normalize){
    r2m = exp(-r2m)
  }else{
    r2m = 1-r2m
  }
  return(r2m)
}

#' @title  Test the significance of lack of fit in SAD/RAD
#'
#' @description  Test the significance of lack of fit in SAD/RAD through Monte-Carlo experiments
#'
#' @param obs the observed abundances as in RAD; no need to be sorted
#' @param model the name of the model as appeared in the name of fit_* function
#' @param params the fitted parameters as returned by fit_* function plus additional arguments like `trunc`
#' @param log logical; whether the goodness/lack of fit is calculated at log-scale; defaults to `TRUE`
#' @param n the number of Monte-Carlo experiments to run for estimating the p-value
#' @param details logical; whether the randomly sampled RADs from Monte-Carlo experiments are returned; defaults to `FALSE`
#' @details Although model selection using AIC can find the best model from a set of SAD models, it is not guaranteed that the best model found is sufficiently good to described the observed SAD/RAD.
#' @details Therefore, `RAD_lack_of_fit_test` compares the fitted SAD with the observed one by calculating the modified R^2 between the observed RAD and the representative RAD given by the fitted SAD.
#' @details The significance of lack of fit is estimated through Monte-Carlo experiments where random RADs are sampled from the fitted SAD and compared to the representative RAD to calculate the baseline empirical distribution for R^2.
#' @details Then the p-value is estimated as the empirical probability of obtaining the observed R^2 or a worse one under the baseline R^2 distribution.
#' @return `RAD_lack_of_fit_test` returns a list of 3 or 4 elements containing the test statistics and other information
#' @return `$R2` the observed R^2
#' @return `$p.value` the p-value of obtaining the observed R^2
#' @return `$representative` the representative RAD given the fitted SAD model
#' @return `$random` optional; the random RADs sampled in the Monte-Carlo experiments
#' @export
#' @rdname RAD_lack_of_fit_test
#'
RAD_lack_of_fit_test = function(obs,model,params,log=TRUE,n=100,details=FALSE,debug=FALSE){
  # variable initialization
  ### future development: argument checkers for meaningful errors and warnings
  numS = length(obs)
  func = paste0("representative_",model)
  # calculate RADs
  predRAD = do.call(func,c(params,list(n=numS,Nrep=n+1,random=FALSE)))
  # calculate R^2
  observed.R2 = r2modified(obs=sort(obs),pred=predRAD[[1]],log=log)
  baseline.R2 = sapply(predRAD[-1],function(x){
    r2modified(obs=x,pred=predRAD[[1]],log=log)
  })
  # estimate p-value
  p.value = sum(observed.R2>=baseline.R2)/n
  # return results
  if(details){
    return(list(R2=observed.R2,p.value=p.value,representative=predRAD[[1]],
                random=predRAD[-1]))
  }else{
    return(list(R2=observed.R2,p.value=p.value,representative=predRAD[[1]]))
  }
}


