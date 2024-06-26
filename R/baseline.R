#' Compute baseline function estimator given the OR parameters
#'
#' Two approaches are used to compute the baseline function estimator:
#'     One uses the weighted (by the odds ratio function) approach
#'     The other uses the iterative approach based on likelihood score
#'
#' @param y a vector (or a matrix) of outcomes.
#' @param x a vector (or a matrix) of covariates.
#' @param parm  a vector (or a matrix) of the OR parameters.
#' @param method "weight" or "iterate" approach for baseline.
#' @param fagg a logical TRUE or FALSE for aggregate the y values
#'               and corresponding baseline function values.
#'
#' @details This function estimates the baseline function
#'                given the OR parameters
#' @return 1. baseline function
#'         2. corresponding y values
#'         3. Location distinct y values in the ordered original y data
#'            if aggregate option is set to TRUE
#'
#' @references Chen, H.Y. (2015). A note ogence of an iterative algorithm for semiparametric odds ratio models.
#'                          Biomatrika, 102, 747-751.
#'
#' @references Chen, H.Y. (2022). Semiparametric Odds Ratio Model and its Application. CRC press.
#'
#'
#' @examples \dontrun{
#' # use the internal data file name dat, a 400x9 data matrix.
#' y=dat[,1:2], x=dat[,3:9];
#' fit=pwlkh(y,x)
#' parm=matrix(fit[[1]],nrow=dim(y)[2])
#' fit1=baseline(y,x,parm)
#' fit2=baseline(y,x,parm,fagg=TRUE)
#' fit3=baseline(y,x,parm,method='iterate')
#' }
#'
#'@export
#'


baseline=function(y,x,parm,method="weight",fagg=TRUE){
  # compute the baseline function estimate with OR parameter estimate
  # y: outcome data
  # x: covariate data
  # parm: regression parameters(if y is a vector, it is lay out column-wise)
  # method: "weight" uses odds ratio as weight to estimate the baseline
  #         "iterate" uses weighted estimate as the initial value,
  #              then iterates by the likelihood score method
  # fagg: whether aggregate the estimate over the distinctive y values
  #

  # 1.initial estimate uses weight
  if(is.vector(y)==TRUE){
    n=length(y)
    p=1
    new=-y*as.vector(matrix(parm,nrow=1)%*%t(x))
  }else{
    n=dim(y)[1]
    p=dim(y)[2]
    new=diag(-y%*%matrix(parm,nrow=dim(y)[2])%*%t(x))  # log(eta(y_i,x_i))
  }
  new=new-max(new)
  logF=new-log(sum(exp(new))) ##keep logarithm of F instead of F to prevent from too close to 0.

  # 2. iterative estimate
  if(method=="iterate"){
    for(i in 1:100){
      #in log-scale
      new=matrix(y,nrow=n)%*%matrix(parm,nrow=p)%*%t(x)  #eta(y_k,x_j)_(nxn). new is a nxn matrix
      logFnew=logF%*%t(rep(1,n))+new                     #log(eta(y,x) dF(y)). logFnew is a nxn matrix
      logFnew=logFnew-rep(1,n)%*%t(apply(logFnew,2,max)) #stablizer is necessary to avoid numeric problem
      logFnew=logFnew-rep(1,n)%*%t(log(apply(exp(logFnew),2,sum)))
                                                         #log(eta(y,x)dF(y)/int eta(y,x)dF(y))
      logFnew=logF-log(apply(exp(logFnew),1,sum))  #-log(int eta(y,x)/int eta(y,x) dF(y) dP_n(x)).
                                                   #logFnew becomes a vector of size n
      logFnew=logFnew-max(logFnew)-log(sum(exp(logFnew-max(logFnew))))

      #print(c(max(logFnew),min(logFnew),sum(exp(logFnew))))

      if(sum(abs(exp(logF)-exp(logFnew)))<1e-5){
        break
      }else{
        #print(sum(abs(exp(logF)-exp(logFnew))))
        logF=logFnew
      }
    }
  }
# Aggregate the weights with the same y observed values.
# The aggregated values will correspond to the candy. So no need to keep y order.
  if(fagg==TRUE){
    ind=c(1:n)
    if(is.vector(y)==TRUE){
      ord=order(y)
      #F=F[ord]
      logF=logF[ord]
      y=y[ord]

      k=1
      for(i in 2:n){
        if(abs(y[i]-y[ind[k]])<1e-8){
          #F[k]=F[k]+F[i]
          logF[k]=logF[k]+log(1+exp(logF[i]-logF[k]))
        }else{
          k=k+1
          #F[k]=F[i]
          logF[k]=logF[i]
          ind[k]=i
        }
      }
    }else{
      ord=do.call(order, as.data.frame(y))
      #F=F[ord]
      logF=logF[ord]
      y=as.matrix(y[ord,])

      k=1
      for(i in 2:n){
        if(sum(abs(y[i,]-y[ind[k],]))<1e-7){
          #F[k]=F[k]+F[i]
          logF[k]=logF[k]+log(1+exp(logF[i]-logF[k]))
        }else{
          k=k+1
          #F[k]=F[i]
          logF[k]=logF[i]
          ind[k]=i
        }
      }

    }

    F=exp(logF) #compute the baseline.
    ###!!! Keep logF in the output is to avoid treating small F as 0, which can be disatrious in this case.
    if(is.vector(y)==TRUE){
      return(list(F[1:k],y[ind[1:k]],logF[1:k],ind[1:k])) # aggregated and ordered
    }else{
      return(list(F[1:k],y[ind[1:k],],logF[1:k],ind[1:k])) # aggregated and ordered
    }
  }else{
    F=exp(logF)
    return(list(F,y,logF)) # nonaggregated and nonordered
  }

}

#fit=pwlkh(x,z)
#parm=matrix(fit[[1]],nrow=dim(x)[2])
#fit1=baseline(x,z,parm,fagg=TRUE)
#fit2=baseline(x,z,parm,method="iterate",fagg=TRUE)

