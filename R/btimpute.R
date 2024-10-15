#' @importFrom graphics par
NULL

#' Implementing bootstrapping multiple simulation-maximization (msm) algorithm by coordinate-wise conditional maximization
#'   1. Create multiple bootstrap data, 
#'   2. run imputation-maximization algorithm for each bootstrap dataset,
#'   3. Run a parallel chain to impute the missing values using the parameter values from
#'        the bootstrap estimates (algorithm).
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
#'  btimpute(dat=bone,miscode=c(-9),nimpute=100)
#' }
#'
#'@export
#'
btimpute=function(dat=dat,miscode=c(-9),method="sp",nem=10,nimpute=10,
               stepsize=0.3,nstep=20,disp=c(3,5)){

  n=dim(dat)[1]
  p=dim(dat)[2]
  misdat=array(1,c(n,p))
  for(k in 1:length(miscode)){
    misdat=misdat-(dat==miscode[k])
  } #assign missing data indicator
  
  impdat=array(0,c(n*nimpute,p))
  #print(dim(impdat))
  bt=sample(x=c(1:n),size=n*nimpute,replace=TRUE)
  for(ni in 1:nimpute){
    impset=(ni-1)*n+c(1:n)
    impdat[impset,]=as.matrix(dat[bt[impset],])
    #print(dim(impdat))
    for(k in 1:p){
      temp=misdat[bt[impset],k]
      obsset=subset(c(1:n),temp==1) # observed data locations
      #print(dim(impdat))
      #print(temp)
      #print(obsset)
      a=table(impdat[(ni-1)*n+obsset,k])
      a_cat=as.numeric(names(a))
      a_freq=as.numeric(a)
      if(sum(1-misdat[bt[impset],k])>0){
        misset=subset(c(1:n),misdat[bt[impset],k]==0) # subset the locations of missing values
        impdat[(ni-1)*n+misset,k]=sample(x=a_cat,prob=a_freq/sum(a_freq),
                                           size=length(misset),replace=TRUE)
      }
    }
  }  # initial imputation for missing values  
  
  #print("here")
  #3. Get the conditional frequencies and impute using conditional draws.
  stheta=array(0,c(nem,p,p-1,nimpute))
  theta=array(0,c(p,p-1,nimpute))
  thetaold=array(0,c(p,p-1,nimpute))
  
  for(iter in 1:nem){
    print(c(iter,nem))
    thetaold=theta
    
    for(k in 1:p){
      print(c(k,k,p))
      for(ni in 1:nimpute){# repeat through imputation
        impset=(ni-1)*n+c(1:n)
        misset=subset(c(1:n),misdat[bt[impset],k]==0) # subset the locations of missing values
        subx=impdat[(ni-1)*n+misset,setdiff(c(1:p),k)] #covariates used for prediction\
        
        if(length(misset)>0){  
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
        }#end of missing
      }#end of multiple imputation
    }#end of k
    

    par(mfrow=disp)
    for(ni in 1:nimpute){
      draw(stheta[,,,ni])
    }
  }#end of em iteration
  
  return(list(theta,impdat,stheta))
}

