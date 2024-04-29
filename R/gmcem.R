#' Implementing Monte Carlo EM algorithm by coordinate-wise conditional maximization
#'
#' @param civ an array civ(n,p,2) with the same number of rows and columns as the
#'          input data and the third dimension has length two to denote lower
#'          and upper bounds of the censoring interval.
#'          For example, if y_k is missing, then lower and upper bounds are -infty
#'          and +infty; We use the observed min and max values to denote them respectively.
#'          If y_k is left censored, the lower bound is set to min and the upper bound
#'          is set to the observed bound. If y_k is right censored, the lower bound
#'          is set to the observed bound and the right bound is set to max.
#'          if y_k is interval censored, the lower bound and upper bound may be
#'          respectively like 4 and 9, meaning the censored value is in the interval (4,9).
#' @param method specific method of estimation, can be
#'          method='sp' for semiparametric likelihood approach (default)
#'          method='pw' for pairwise likelihood approach
#'          method='pm' for permutation likelihood approach
#' @param nem the number of EM steps
#' @param nimpute the number of Monte Carlo copies for each missing value
#' @param stepsize controls the step size of each M-step update.
#' @param nstep controls the number of steps before the change of step size
#'
#' @details This function maximizes the conditional likelihood coordinate-wise
#'          by the Monte Carlo EM algorithm. In M-step, all the OR parameters are
#'          estimated by the semiparametric OR model using the semiparametric likelihood
#'          approach based on the current imputed data. In MC-step, baseline is first
#'          obtained. Imputing missing values are done coordinate-wise by
#'          the sampling approach
#'
#' @return  1. parameter estimate
#'          2. one copy of the imputed data
#'
#' @references Chen, H.Y. (2022). Semiparametric Odds Ratio Model and its Application. CRC press.
#'
#'
#' @examples \dontrun{
#'  gmcem(civ=extbone, nem=100, stepsize=0.6,nstep=5)
#' }
#'
#'@export
#'

gmcem=function(civ,method="sp",nem=10,nimpute=1,stepsize=0.5,nstep=20){
  #1. preparation
  #1.1. Find the missing/censoring data indicators
  n=dim(civ)[1]
  p=dim(civ)[2]
  censor=1.0*(civ[,,1]==civ[,,2]) #missing/censoring indicator

  #1.2. Set the imputed data set
  impdat=array(0,c(n*nimpute,p+1))
  for(ni in 1:nimpute){
    dstart=(ni-1)*n+1
    dend=ni*n
    impdat[dstart:dend,1:p]=as.matrix(civ[,,1])
    impdat[dstart:dend,p+1]=ni  # add a variable index the datasets
  }
  #1.3. Transform the original data into locations.
  #     The original observed data values are retained in the impdat.
  for(k in 1:p){
    if(sum(1-censor[,k])>0){
      obsset=subset(c(1:n),censor[,k]==1) # observed data locations
      a=table(civ[obsset,k,1])
      a_cat=as.numeric(names(a)) # categories of the observed data
      a_freq=as.numeric(a)       # frequency of each category
      censet=setdiff(c(1:n),obsset) # censored data locations in the sample
 
      civ[censet,k,1]=insertloc(a_cat,civ[censet,k,1],method="nearestleft")[[1]]
              # Convert lower bound values of censoring intervals into 
              # locations in the ordered completely observed data values (a_cat)
              # This conversion makes it easy to sample in the imputation step!!
      civ[censet,k,2]=insertloc(a_cat,civ[censet,k,2],method="nearestright")[[1]] 
              # Convert upper bound values of censoring intervals into 
              # locations in the ordered completely observed data values (a_cat)
              # This conversion makes it easy to sample in the imputation step!!

      for(j in 1:length(censet)){
        tx=a_cat[civ[censet[j],k,1]:civ[censet[j],k,2]]
        tprob=a_freq[civ[censet[j],k,1]:civ[censet[j],k,2]]
        tprob=tprob/sum(tprob)
        impdat[c(0:(nimpute-1))*n+censet[j],k]=sample(x=tx,prob=tprob,size=nimpute,replace=TRUE)
      }
    }  
  }

  #3. Get the conditional frequencies and impute using random draws.
  stheta=array(0,c(nem,p,p-1))
  theta=array(0,c(p,p-1))
  for(iter in 1:nem){
    print(c(iter,nem))
    for(k in 1:p){

      print(c(k,k,p))
      thetaold=theta[k,]
#      if(sum(1-misdat[,k])>0){

        #a. M-step
        #print('M-step')
        if(method=='pw'){
          fit=pwlkh(impdat[,k],impdat[,setdiff(c(1:p),k)])
        }else{
          if(method=='sp'){
            fit=splkh(impdat[,k],impdat[,setdiff(c(1:p),k)])
          }else{
            fit=pmlkh(cbind(impdat[,k],impdat[,setdiff(c(1:p),k)]),c(1,p-1))
          }
        }
        theta[k,]=fit[[1]]

        if(iter>nstep){stepsize=1}
        theta[k,]=thetaold+stepsize*(theta[k,]-thetaold)

    print(theta[k,])
    print(sum(abs(theta[k,]-thetaold)))

        # b. MC-step
        if(sum(1-censor[,k])>0){
          #print('MC-step')
          base=baseline(impdat[,k],impdat[,setdiff(c(1:p),k)],parm=theta[k,],method="iterate",fagg=TRUE)
          #base=baseline(impdat[,k],impdat[,setdiff(c(1:p),k)],parm=theta[k,],fagg=TRUE)

          censet=subset(c(1:n),censor[,k]==0) # subset the locations of missing values
          subx=impdat[censet,setdiff(c(1:p),k)]
          pred=cprob(y=base[[2]],x=subx,parm=theta[k,],logF=base[[3]])
          
         # print(base[[2]])
          for(j in 1:length(censet)){
            ysubset=c(civ[censet[j],k,1]:civ[censet[j],k,2])  
           # print(c(j,censet[j],censet[j]))
          #  print(ysubset)
            
            predy=base[[2]][ysubset]         # This is also the same as a_cat
            
            logprob=pred[[2]][ysubset,j]     # this is the log-probabilities
            logprob=logprob-rep(max(logprob),length(ysubset))  #To rescale to make it close to 0.
                                                               #This is necessary to avoid numeric problem.
            logprob=logprob-rep(log(sum(exp(logprob))),length(ysubset))
            predprob=exp(logprob)
 
        #  if(length(predy)>1){
            impdat[c(0:(nimpute-1))*n+censet[j],k]=sample(x=predy,prob=predprob,size=nimpute,replace=TRUE)
        #  }else{
        #    impdat[c(0:(nimpute-1))*n+censet[j],k]=rep(predy,nimpute) #if interval includes only one value, just repeat
        #  }
          }
       }
    }
    stheta[iter,,]=theta
    draw(stheta)

  }

  return(list(theta,impdat,stheta))
}

