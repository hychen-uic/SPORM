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

cdimpute=function(dat=dat,miscode=c(-9),method="pw",niter=10){
  #1. Find the missing data indicators
  n=dim(dat)[1]
  p=dim(dat)[2]
  misdat=array(1,c(n,p))
  for(k in 1:length(miscode)){
    misdat=misdat-(dat==miscode[k])
  } # missing data indicators

  #2. Get the marginal frequencies and impute using random draws.
  impdat=dat
  for(k in 1:p){
    a=table(dat[,k])
    a_cat=as.numeric(names(a))
    a_freq=as.numeric(a)
    if(sum(1-misdat[,k])>0){
      len=length(a_cat)
      misset=order(misdat[,k])[1:a_freq[1]] # subset the locations of missing values
      imp=sample(x=a_cat[2:len],prob=a_freq[2:len]/sum(a_freq[2:len]),size=a_freq[1],replace=TRUE)
      impdat[misset,k]=imp # fill in the imputed values by random draws
    }
  }

  #3. Get the conditional frequencies and impute using random draws.
  impdat=as.matrix(impdat)
  for(iter in 1:niter){
    print(c(iter,iter,niter))
    for(k in 1:p){
      print(c(k,k,p))
      a=table(dat[,k])
      a_cat=as.numeric(names(a))
      a_freq=as.numeric(a)
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
        print(fit[[1]])
        base=baseline(impdat[,k],impdat[,setdiff(c(1:p),k)],parm=fit[[1]],fagg=TRUE)
        print(base)
        len=length(a_cat)
        misset=order(misdat[,k])[1:a_freq[1]] # subset the locations of missing values
        subx=impdat[misset,setdiff(c(1:p),k)]
        pred=cprob(y=base[[2]],x=subx,parm=fit[[1]],F=base[[1]])
        #print(dim(pred[[1]]));print(dim(subx))
        imp=rep(0,length(misset))
        for(j in 1:length(misset)){
          imp[j]=sample(x=base[[2]],prob=pred[[1]][,j],size=1)
        }
        impdat[misset,k]=imp
        #        print(table(dat[,k]))
        #        print(table(impdat[,k]))
      }
    }
  }
  return(list(impdat))
}

