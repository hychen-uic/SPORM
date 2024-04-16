#' Compute conditional probability for a given covariate x, outcome y, OR parameter theta, and baseline function F.
#'
#' @param y a vector (or a matrix) of all possible outcome values.
#' @param x a vector of covariates.
#' @param parm  a vector (or a matrix) of the OR parameters.
#' @param logF logarithm (to keep from treating small probability as 0) of baseline probabilities for all possible y values
#'
#' @details This function estimates the predictive probabilities for
#'               given parameters theta and baseline function
#' @return 1. all the predictive probabilities for a given covariate x.
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
#' theta=matrix(fit[[1]],nrow=dim(y)[2])
#' fit=baseline(y,x,parm,method='iterate')
#' fit1=cprob(y,x,parm=theta,F=fit[[1]])
#' predprob=fit1[[1]]
#' }
#'
#'@export
#'

cprob=function(y,x,parm,logF){

  # compute the all the predictive probabilities for a given x value
  # y: all possible outcome values
  # x: given covariate values.
  # parm: regression parameters(if y is a vector, it is lay out column-wise)
  # logF: logarithm of the baseline function. Use log-scale to avoid treating very small F as 0.
  # logpred: logarithm of the predictive probabilities
  #
  if(is.vector(y)==TRUE){
    n=length(y)
    p=1
  }else{
    n=dim(y)[1]
    p=dim(y)[2]
  }
  if(is.vector(x)==TRUE){
    nmiss=length(x)
  }else{
    nmiss=dim(x)[1]
  }
  new=matrix(y,nrow=n)%*%matrix(parm,nrow=p)%*%t(x)  #log(eta(y_k,x_j)_(nxnmiss)). x_(nmissxp)
  logpred=logF%*%t(rep(1,nmiss))+new           # log(eta(y,x) dF(y)). logpred is a nxnmiss matrix
                                               # This is in log-scale to avoid treating very small pred as 0.
  logpred=logpred-rep(1,n)%*%t(apply(logpred,2,max))  #To rescale to make it close to 0.
                                                      #This is necessary to avoid numeric problem.
  logpred=logpred-rep(1,n)%*%t(log(apply(exp(logpred),2,sum)))
                        # log(eta(y,x) dF(y)/int eta(y,x) dF(y)). pred is a nxnmiss matrix
  pred=exp(logpred)     # the predictive probabilities for each missing value

  return(list(pred))
}


