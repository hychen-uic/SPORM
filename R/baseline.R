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

  if(is.vector(y)==TRUE){
    n=length(y)
    F=rep(0,n)
    # initial estimate uses weight
    F=exp(-y*as.vector(x%*%parm))
    F=F/sum(F)
  }else{
    # initial estimate uses weight
    n=dim(y)[1]
    F=diag(exp(y%*%parm%*%t(x)))
    #F=rep(0,n)
    #for(i in 1:n){
    #  F[i]=exp(-as.numeric(y[i,]%*%parm%*%x[i,]))
    #  }
    F=F/sum(F)
  }

  if(method=="iterate"){
    for(i in 1:50){
      if(is.vector(parm)==TRUE){
        new=exp(y%*%t(parm)%*%t(x))
      }else{
        new=exp(y%*%parm%*%t(x))
      }
      Fnew=as.vector(F%*%new)
      Fnew=Fnew/diag(new)
      Fnew=Fnew/sum(Fnew)
      if(sum(abs(F-Fnew))<1e-6){
        #print(c(i,i,i))
        break
      }else{
        F=Fnew
      }
    }
  }
# Aggregate the weights with the same y observed values.
# The aggregated values will correspond to the candy. So no need to keep y order.
  if(fagg==TRUE){
    ind=c(1:n)
    if(is.vector(y)==TRUE){
      ord=order(y)
      F=F[ord]
      y=y[ord]

      k=1
      for(i in 2:n){
        if(abs(y[i]-y[ind[k]])<1e-8){
          F[k]=F[k]+F[i]
        }else{
          k=k+1
          F[k]=F[i]
          ind[k]=i
        }
      }
    }else{
      ord=do.call(order, as.data.frame(y))
      F=F[ord]
      y=y[ord,]

      k=1
      for(i in 2:n){
        if(sum(abs(y[i,]-y[ind[k],]))<1e-7){
          F[k]=F[k]+F[i]
        }else{
          k=k+1
          F[k]=F[i]
          ind[k]=i
        }
      }

    }

    if(is.vector(y)==TRUE){
      return(list( F[1:k],y[ind[1:k]],ind[1:k])) # aggregated and ordered
    }else{
      return(list(F[1:k],y[ind[1:k],],ind[1:k])) # aggregated and ordered
    }
  }else{
    return(list(F,y)) # nonaggregated and nonordered
  }

}

#fit=pwlkh(x,z)
#parm=matrix(fit[[1]],nrow=dim(x)[2])
#fit1=baseline(x,z,parm,fagg=TRUE)
#fit2=baseline(x,z,parm,method="iterate",fagg=TRUE)

