# microSAD: Sampling-explicit fitting of microbial SAD
`microSAD` is an open-source R package designed for testing Species Abundance Distribution (SAD) models within microbial communities. It addresses unique challenges inherent in fitting microbial SADs, including sampling processes in the 16S rRNA sequencing pipeline, variability in 16S rRNA gene copy numbers among species, large species abundance, and small probability mass that may affect model fitting accuracy in commonly used R packages like `sads`.

Specifically, `microSAD` has useful features related to microbial SAD data:

1. explicit Poisson/negative-binomial sampling error incorporated to address general sampling effort and bias

2. enable OTU-specific sampling effort/bias correction (e.g., 16S rRNA GCN variation among OTUs)

3. convenient functions for fitting SAD, predicting RAD, and testing the significance of lack of fit with good precision and speed

At present, `microSAD` offers the capability to test four base SAD models: logseries, lognormal, power law, and powerbend, along with corresponding compound SAD models featuring Poisson or negative binomial error structures.

While primarily intended for microbial SAD data analysis, `microSAD` is also applicable for testing SAD models in animal and plant communities.

Detailed description of the methods and analyses are available in the preprint: `A unifying model of species abundance distribution`

doi: --

## System requirements
The package has been tested on using R version 4.3.2 (2023-10-31).
The following R packages are required:
`sads`,`pracma`

## Installation
`microSAD` can be installed using the following command in R
```
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github(repo = "wu-lab-uva/microSAD")
```
If not installed, the required packages can be installed from CRAN
```
list.of.packages <- c("pracma", "sads")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages,repos = "http://cran.us.r-project.org", INSTALL_opts = '--no-lock')
```
## Data format
Similar to `sads`, `microSAD` fits SAD from a numeric vector of species abundances (e.g., number of sequence reads for an OTU). 
Instead of saving the fitted results in an `sads` object, the fitting functions in `microSAD` saves the fitted model in an `mle2` object.

## Demo
Once installed, a small demo can be run in R to check if `microSAD` operates properly:
```
### load demo data
library(microSAD)
data2fit = readRDS(system.file("extdata/Demo","demo_SAD.RDS",package = "microSAD",mustWork = TRUE))
sprintf("Number of sequence reads in the SAD data: %d",sum(data2fit))
sprintf("Number of species in the SAD data: %d",length(data2fit))

### fit Poisson lognormal
start = Sys.time()
res = list()
try({
    this.res = fit_poilog(x=data2fit,trunc=0,control=list(maxit=2000))
    res$poilog = list(family="lnorm",model="lnorm",sampling="Poisson",full.model= "poilog",
                      coef=this.res@coef,fullcoef=this.res@fullcoef,AIC=AIC(this.res),
                      vcov = this.res@vcov,convergence=this.res@details$convergence)
})
### fit logseries
try({
    this.res = fitls(x=data2fit,trunc=0)
    res$ls= list(family="ls",model="ls",sampling="none",full.model= "ls",
                 coef=this.res@coef,fullcoef=this.res@fullcoef,AIC=AIC(this.res),
                 vcov = this.res@vcov,convergence=this.res@details$convergence)
})
# if the previous attempt fails, fit logseries again using a starting value of 1 for the parameter 
if(is.null(res$ls)){
    try({
      this.res = fitls(x=data2fit,trunc=0,start.value=1)
      res$ls= list(family="ls",model="ls",sampling="none",full.model= "ls",
                   coef=this.res@coef,fullcoef=this.res@fullcoef,AIC=AIC(this.res),
                   vcov = this.res@vcov,convergence=this.res@details$convergence)
    })
}
### fit Poisson logeries
try({
    this.res = fit_poipowbend(x=data2fit,trunc=0,fixed=list(s=1),control=list(maxit=2000))
    res$poils= list(family="ls",model="ls",sampling="Poisson",full.model= "poils",
                    coef=this.res@coef,fullcoef=this.res@fullcoef,AIC=AIC(this.res),
                    vcov = this.res@vcov,convergence=this.res@details$convergence)
})
### fit powerbend
# fit powerbend with BFGS
try({
    this.res = fitpowbend(x=data2fit,trunc=0,control=list(maxit=2000))
    res$powbend = list(family="powbend",model="powbend",sampling="none",full.model= "powbend",
                       coef=this.res@coef,fullcoef=this.res@fullcoef,AIC=AIC(this.res),
                       vcov = this.res@vcov,convergence=this.res@details$convergence)
},silent = TRUE)
# fit powerbend with Nelder-Mead if BFGS fails
if(is.null(res$powbend)){
    try({
      this.res = fitpowbend(x=data2fit,trunc=0,start.value = c(1,7),method="Nelder-Mead",control=list(maxit=2000))
      res$powbend = list(family="powbend",model="powbend",sampling="none",full.model= "powbend",
                         coef=this.res@coef,fullcoef=this.res@fullcoef,AIC=AIC(this.res),
                         vcov = this.res@vcov,convergence=this.res@details$convergence)
    })
}

# if the previous attempt fails, fit powerbend using logseries model parameters as starting values if they exist 
if(is.null(res$powbend)){
    try({
      this.res = fitpowbend(x=data2fit,trunc=0,start.value = c(1,-log(-log(res$ls$fullcoef["N"]/sum(res$ls$fullcoef)))),method="Nelder-Mead",control=list(maxit=2000))
      res$powbend = list(family="powbend",model="powbend",sampling="none",full.model= "powbend",
                         coef=this.res@coef,fullcoef=this.res@fullcoef,AIC=AIC(this.res),
                         vcov = this.res@vcov,convergence=this.res@details$convergence)
    })
}
#### fit Poisson powerbend
#fit Poisson powerbend using powerbend model parameters as starting values if they exist
if(is.null(res$powbend)){
    start.value = list(s=1,minlogLambda=7,eta=0.001)
}else{
    start.value = list(s=unname(res$powbend$fullcoef["s"]),
                       minlogLambda=unname(res$powbend$fullcoef["oM"]),eta=0.001)
}
try({
    this.res = fit_poipowbend(x=data2fit,trunc=0,start.value = start.value,control=list(maxit=2000))
    res$poipowbend = list(family="powbend",model="powbend",sampling="Poisson",full.model= "poipowbend",
                          coef=this.res@coef,fullcoef=this.res@fullcoef,AIC=AIC(this.res),
                          vcov = this.res@vcov,convergence=this.res@details$convergence)
})

# fit Poisson powerbend using Poisson logseries model parameters as starting values if they exist
if(is.null(res$poils)){
    start.value = list(s=1,minlogLambda=7,eta=0.001)
}else{
    start.value = as.list(res$poils$fullcoef)
}
try({
    this.res = fit_poipowbend(x=data2fit,trunc=0,start.value = start.value,control=list(maxit=2000))
    this.overwrite=is.null(res$poipowbend)
    if(!this.overwrite){
      this.overwrite = (AIC(this.res)<res$poipowbend$AIC)
    }
    if(this.overwrite){
      res$poipowbend = list(family="powbend",model="powbend",sampling="Poisson",full.model= "poipowbend",
                            coef=this.res@coef,fullcoef=this.res@fullcoef,AIC=AIC(this.res),
                            vcov = this.res@vcov,convergence=this.res@details$convergence)
    }
})
### fit power law
try({
      this.res = fitpower(x=data2fit,trunc=0,control=list(maxit=2000))
      res$power = list(family="power",model="power",sampling="None",full.model= "power",
                       coef=this.res@coef,fullcoef=this.res@fullcoef,AIC=AIC(this.res),
                       vcov = this.res@vcov,convergence=this.res@details$convergence)
})
### fit Poisson power law 
# fit Poisson power law using power law model parameters as starting values if they exist
if(is.null(res$power$fullcoef['s'])){
      start.value = list()
}else { 
      start.value = list(s=res$power$fullcoef['s'],eta=1)
}
try({
    this.res = fit_poipower(x=data2fit,trunc=0,start.value = start.value, control=list(maxit=2000))
    res$poipower = list(family="power",model="power",sampling="Poisson",full.model= "poipower",
                      coef=this.res@coef,fullcoef=this.res@fullcoef,AIC=AIC(this.res),
                      vcov = this.res@vcov,convergence=this.res@details$convergence)
})

#show AIC and fitted model parameters
cat("AIC\n")
print(sapply(res,function(x){x$AIC}))
cat("Parameters\n")
print(sapply(res,function(x){x$fullcoef}))

end = Sys.time()
sprintf("It took %f minutes to complete the model fitting.", difftime(end,start,units = "mins"))
```
If everything works well, the following output should be expected in 5 to 10 minutes. There may be small differences in the fitted parameters and AIC due to stochastic procedures in the fitting algorithm but they outcomes should be highly similar, if not identical. The time required for fitting may vary depending on your hardware settings.
```
[1] "Number of sequence reads in the SAD data: 23610"
[1] "Number of species in the SAD data: 1741"

AIC
    poilog         ls      poils    powbend poipowbend      power   poipower 
  6262.795   8006.105   8008.105   6330.000   6259.625   6328.019   6263.385 
Parameters
$poilog
       mu       sig 
-85.34373  11.51679 

$ls
       N    alpha 
23610.00   433.56 

$poils
           s minlogLambda          eta 
1.000000e+00 2.318513e+01 4.642665e-09 

$powbend
        s        oM 
 1.906497 12.775317 

$poipowbend
           s minlogLambda          eta 
1.632911e+00 4.006902e+01 1.653405e-14 

$power
       s 
1.906799 

$poipower
           s          eta 
1.6494890330 0.0002023111 

[1] "It took 4.494885 minutes to complete the model fitting."
```