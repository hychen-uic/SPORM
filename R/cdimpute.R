#' Impute missing values by random draw from the conditional distributions
#'
#' This is a repeated imputation by the conditional distribution
#'
#' @param dat a data matrix with possible missing values
#' @param miscode a set of values used for denoting missing values
#' @param method specific method of estimation, can be
#'        method='pw' for pairwise likelihood approach
#'        method='sp' for semiparametric likelihood approach
#'        method='pm' for permutation likelihood approach
#' @param niter the number of iterations to produce the final imputed data set
#' @param nimpute the number of imputed copies for each missing value
#'
#' @details This function imputes missing data by random draw from
#' the conditional distribution of the variable conditional on other variables.
#' in each iteration, the conditional model is estimated by
#'   the semiparametric OR model using the semiparametric likelihood approach
#'   based on the current imputed data
#'
#' @return  1. One imputed data set
#'
#' @references Chen, H.Y. (2022). Semiparametric Odds Ratio Model and its Application. CRC press.
#'
#'
#' @examples \dontrun{
#'  cdimpute(dat=bone,miss=c(-9))
#' }
#'
#'@export
#'

cdimpute=function(dat=dat,miscode=c(-9),method="sp",niter=10,nimpute=5){
  #1. Find the missing data indicators
  n=dim(dat)[1]
  p=dim(dat)[2]
  misdat=array(1,c(n,p))
  for(k in 1:length(miscode)){
    misdat=misdat-(dat==miscode[k])
  } # missing data indicators

  #2. Get the conditional frequencies and impute using random draws.
  impdat=array(0,c(n*nimpute,p+1))
  for(ni in 1:nimpute){
    dstart=(ni-1)*n+1
    dend=ni*n
    impdat[dstart:dend,1:p]=as.matrix(dat)
    impdat[dstart:dend,p+1]=ni  # add a variable index the datasets
  }
  for(k in 1:p){
    obsset=subset(c(1:n),misdat[,k]==1) # observed data locations
    a=table(dat[obsset,k])
    a_cat=as.numeric(names(a))
    a_freq=as.numeric(a)
    if(sum(1-misdat[,k])>0){
      misset=subset(c(1:n),misdat[,k]==0) # subset the locations of missing values
      imp=sample(x=a_cat,prob=a_freq/sum(a_freq),size=length(misset)*nimpute,replace=TRUE)
      imp=matrix(imp, ncol=nimpute)
      for(ni in 1:nimpute){
        dstart=(ni-1)*n+1
        dend=ni*n
        impdat[(ni-1)*n+misset,k]=imp[,ni] # fill in the imputed values by random draws
      }
    }
  }

  #3. Get the conditional frequencies and impute using random draws.
  for(iter in 1:niter){
    print(c(iter,iter,niter))
    for(k in 1:p){
      print(c(k,k,p))
      if(sum(1-misdat[,k])>0){
        if(method=='pw'){
          fit=pwlkh(impdat[,k],impdat[,setdiff(c(1:p),k)])
        }else{
          if(method=='sp'){
            fit=splkh(impdat[,k],impdat[,setdiff(c(1:p),k)])
          }else{
            fit=pmlkh(cbind(impdat[,k],impdat[,setdiff(c(1:p),k)]),c(1,p-1))
          }
        }
        base=baseline(impdat[,k],impdat[,setdiff(c(1:p),k)],parm=fit[[1]],fagg=TRUE)

        misset=subset(c(1:n),misdat[,k]==0) # subset the locations of missing values
        subx=impdat[misset,setdiff(c(1:p),k)]
        pred=cprob(y=base[[2]],x=subx,parm=fit[[1]],F=base[[1]])

        imp=array(0,c(length(misset),nimpute))
        for(j in 1:length(misset)){
          imp[j,]=sample(x=base[[2]],prob=pred[[1]][,j],size=nimpute,replace=TRUE)
        }
        for(ni in 1:nimpute){
          dstart=(ni-1)*n+1
          dend=ni*n
          impdat[(ni-1)*n+misset,k]=imp[,ni]          }
        }
       if(k==1){
         theta=fit[[1]]
       }else{if(k<p){
         theta=c(theta,fit[[1]][k:(p-1)])
         }}
     }
  }
  print(c(min(theta),max(theta)))
  return(list(impdat,theta))
}

#fit=cdimpute(bone, miscode=c(-9),niter=20,nimpute=200)
