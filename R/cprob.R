#' Compute conditional probability for a given covariate x, outcome y, OR parameter theta, and baseline function F.
#'
#' @param y a vector (or a matrix) of all possible outcome values.
#' @param x a vector of covariates.
#' @param parm  a vector (or a matrix) of the OR parameters.
#' @param F baseline probabilities for all possible y values
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

cprob=function(y,x,parm,F){
  # compute the all the predictive probabilities for a given x value
  # y: all possible outcome values
  # x: given covariate value
  # parm: regression parameters(if y is a vector, it is lay out column-wise)
  # F: baseline function
  # pred: predictive probabilities corresponding to all possible y values.
  #
  if(is.vector(y)==TRUE){
    top=y%*%matrix(parm,nrow=1)%*%t(x)
  }else{
    top=y%*%matrix(parm,nrow=dim(y)[2])%*%t(x)
  }
  pred=diag(log(F))%*%top
#  print(pred)
  pred=pred-rep(1,dim(pred)[1])%*%t(apply(pred,2,max))
      # subtract different rows by different constants
#  print(pred)
  pred=pred%*%diag(1/apply(pred,2,sum))
#  print('here')
#  print(parm)
#  print(F)
#  print(pred)

  return(list(pred))
}


