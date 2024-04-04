#' Implementing Monte Carlo EM algorithm by coordinate-wise conditional maximization
#'
#' @param dat a data matrix with possible missing values
#' @param miscode a set of values used for denoting missing values
#' @param method specific method of estimation, can be
#'         method='sp' for semiparametric likelihood approach (default)
#'        method='pw' for pairwise likelihood approach
#'        method='pm' for permutation likelihood approach
#' @param nem the number of EM steps
#' @param nimpute the number of Monte Carlo copies for each missing value
#' @param nseq the number of repeated sampling steps in each MC-step
#'
#' @details This function maximizes the conditional likelihood coordinate-wise
#' by the Monte Carlo EM algorithm. In M-step, all the OR parameters are
#' estimated by the semiparametric OR model using the semiparametric likelihood
#' approach based on the current imputed data. In MC-step, baseline is first
#' obtained. Imputing missing values are done coordinate-wise by
#' the sampling approach
#'
#' @return  1. parameter estimate
#'          2. one copy of the imputed data
#'
#' @references Chen, H.Y. (2022). Semiparametric Odds Ratio Model and its Application. CRC press.
#'
#'
#' @examples \dontrun{
#'  mcem(dat=bone,miss=c(-9),nimpute=100)
#' }
#'
#'@export
#'

mcem=function(dat=dat,miscode=c(-9),method="sp",nem=10,nimpute=5,nseq=10){
  #1. Find the missing data indicators
  n=dim(dat)[1]
  p=dim(dat)[2]
  misdat=array(1,c(n,p))
  for(k in 1:length(miscode)){
    misdat=misdat-(dat==miscode[k])
    } # missing data indicators

  #2. Get the marginal frequencies and impute using random draws.
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
  theta=array(0,c(p,p-1))
  for(iter in 1:nem){
    print(c(iter,nem))
    for(k in 1:p){
      #a. M-step
      #print('M-step')
      print(c(k,k,p))
      thetaold=theta[k,]
      if(sum(1-misdat[,k])>0){
        if(method=='pw'){
          fit=pwlkh(impdat[,k],impdat[,setdiff(c(1:p),k)])
        }else{
          if(method=='sp'){
            fit=splkh(impdat[,k],impdat[,setdiff(c(1:p),k)])
          }else{
            fit=pmlkh(cbind(impdat[,k],impdat[,setdiff(c(1:p),k)]),c(1,p-1))
          }
        theta[k,]=fit[[1]]
        }

   #    }
    #print(sum(abs(thetaold-theta[k,])))
    #print(c(iter,iter,iter,nem))
    #print(c(min(theta[k,]),max(theta[k,]),sum(abs(theta[k,]-thetaold))))
    print(theta[k,])
    print(sum(abs(theta[k,]-thetaold)))
    #thetaold=theta

    #print('MC-step')
    #for(nrep in 1:nseq){
    #  print(c(nrep,nseq))
    #for(k in 1:p){
    #  print(c(k,k,p))
    #  if(sum(1-misdat[,k])>0){
        base=baseline(impdat[,k],impdat[,setdiff(c(1:p),k)],parm=theta[k,],method="iterate",fagg=TRUE)

        #plot(base[[2]],base[[1]])

        misset=subset(c(1:n),misdat[,k]==0) # subset the locations of missing values
        subx=impdat[misset,setdiff(c(1:p),k)]

        pred=cprob(y=base[[2]],x=subx,parm=theta[k,],F=base[[1]])
        #print(pred[[1]])

        imp=array(0,c(length(misset),nimpute))
        for(j in 1:length(misset)){
          imp[j,]=sample(x=base[[2]],prob=pred[[1]][,j],size=nimpute,replace=TRUE)
        }
        for(ni in 1:nimpute){
          dstart=(ni-1)*n+1
          dend=ni*n
          impdat[(ni-1)*n+misset,k]=imp[,ni]
        }
      }
    }
  }

  return(list(theta,impdat))
}

