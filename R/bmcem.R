#' @importFrom graphics par
NULL

#' Implementing multiple simulation-maximization (msm) algorithm by coordinate-wise conditional maximization
#' Do not average of parameters in each iteration. Each imputation and maximization step is independent 
#' of the others.
#'
#' @param dat a data matrix with possible missing values
#' @param miscode a set of values used for denoting missing values
#' @param method specific method of estimation, can be
#'         method='sp' for semiparametric likelihood approach (default)
#'         method='pw' for pairwise likelihood approach
#'         method='pm' for permutation likelihood approach
#' @param nem the number of EM steps
#' @param nimpute the number of Monte Carlo copies for each missing value
#' @param stepsize controls the step size of each M-step update.
#' @param nstep controls the number of steps before the change of step size
#' @param disp controls the layout of the parameter history plot: a two-dimensional
#'        vector of integers in par, mfrow=disp=c(3,5)
#'        
#' @details This function maximizes the conditional likelihood coordinate-wise
#' by the Monte Carlo EM algorithm. In M-step, all the OR parameters are
#' estimated by the semiparametric OR model using the semiparametric likelihood
#' approach based on the current imputed data. In MC-step, baseline is first
#' obtained. Imputing missing values are done coordinate-wise by
#' the sampling approach
#'
#' @return  1. parameter estimate for theta
#'          2. one copy of the imputed data
#'          3. parameter estimates for theta in successful iterations.
#'
#' @references Chen, H.Y. (2022). Semiparametric Odds Ratio Model and its Application. CRC press.
#'
#'
#' @examples \dontrun{
#'  amcem(dat=bone,miss=c(-9),nimpute=100)
#' }
#'
#'@export
#'
bmcem=function(dat=dat,miscode=c(-9),method="sp",nem=10,nimpute=10,
               stepsize=0.3,nstep=20,disp=c(3,5)){

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
      for(ni in 1:nimpute){
        impdat[(ni-1)*n+misset,k]=sample(x=a_cat,prob=a_freq/sum(a_freq),
                                      size=length(misset),replace=TRUE)
      }
    }
  }

  #3. Get the conditional frequencies and impute using conditional draws.
  stheta=array(0,c(nem,p,p-1,nimpute))
  theta=array(0,c(p,p-1,nimpute))
  thetaold=array(0,c(p,p-1,nimpute))
  
  for(iter in 1:nem){
    print(c(iter,nem))
    thetaold=theta
    
    for(k in 1:p){
      print(c(k,k,p))
      
      if(sum(1-misdat[,k])>0){ 
  
        misset=subset(c(1:n),misdat[,k]==0) # subset the locations of missing values
        for(ni in 1:nimpute){# repeat through imputation
          impset=(ni-1)*n+c(1:n)
          subx=impdat[(ni-1)*n+misset,setdiff(c(1:p),k)] #covariates used for prediction
          
      #a. M-step
        
          if(method=='pw'){
            fit=pwlkh(impdat[impset,k],impdat[impset,setdiff(c(1:p),k)])
          }else{
            if(method=='sp'){
              fit=splkh(impdat[impset,k],impdat[impset,setdiff(c(1:p),k)])
            }else{
              fit=pmlkh(cbind(impdat[impset,k],impdat[impset,setdiff(c(1:p),k)]),c(1,p-1))
            }
          }
        
          theta[k,,ni]=fit[[1]]
          if(iter>nstep){stepsize=1}
          theta[k,,ni]=thetaold[k,,ni]+stepsize*(theta[k,,ni]-thetaold[k,,ni])
          stheta[iter,k,,ni]=theta[k,,ni] #record the iteration history for theta
          
          #print(theta[k,,ni])
          #print(sum(abs(theta[k,,ni]-thetaold[k,,ni])))
          
          # b. MC-step
          #print('MC-step')
        
          base=baseline(impdat[impset,k],impdat[impset,setdiff(c(1:p),k)],parm=theta[k,,ni],method="iterate",fagg=TRUE)
          pred=cprob(y=base[[2]],x=subx,parm=theta[k,,ni],logF=base[[3]])
     
          for(j in 1:length(misset)){
            impdat[(ni-1)*n+misset[j],k]=sample(x=base[[2]],prob=pred[[1]][,j],size=1,replace=TRUE)
          }  
        }
      }#end of missing
    }#end of k
    
    par(mfrow=disp)
    for(ni in 1:nimpute){
      draw(stheta[,,,ni])
    }
  }#end of em iteration
  
  return(list(theta,impdat,stheta))
}

