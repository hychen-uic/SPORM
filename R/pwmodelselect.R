#'
#'  Penalized model selection by L^1 penalized pairwise pseudo-likelihood
#'   with AIC type of criterion for tuning parameter selection
#'
#' This approach selects an odds ratio model
#'    using the penalized pairwise pseudo-likelihood approach
#'    with AIC for tuning parameter selection.
#'
#' @param y a vector of outcomes
#' @param w a vector of covariates
#' @param lamb a vector of penalty values for model selection,
#'               each penalty value determine a model.
#'
#' @details The final selected model is determined by refitting
#'             using AIC criterion.
#' @return a set of selected parameters to determine the model, AIC values,
#'           and the structure of the model
#'
#' @references Chen, H.Y. (2022). Semiparametric Odds Ratio Model and its Application. CRC press.
#'
#' @references Liang, K.Y. and Qin, J. (2000). Regression analysis under non-standard situations:
#'             a pairwise pseudolikelihood approach, Journal of the Royal Statistical Society, Ser. B,
#'             62, 773-786.
#'
#' @examples \dontrun{
#' # use the internal data file name dat, a 400x9 data matrix.
#' y=dat[,1:2], x=dat[,3:9];lambda=c(10,50,200,500,800)
#' pwmodelselect(y,w,lamb)
#' }
#'
#'@export
#'

pwmodelselect=function(y,w,lamb){

    loglkh=rep(0,length(lamb))
    AIC=rep(0,length(lamb))

    if(is.vector(y)==TRUE){
      n=length(y)
    }else{
      n=dim(y)[1]
    }

    fit=pwhdlkh(y,w,lambda=lamb) #in SPORM
    for(k in 1:length(lamb)){
      loglkh[k] <- pwlkhfix(y,w,fixstruct=fit[[2]][,,k])[[3]] #in SPORM
      AIC[k]=-loglkh[k]/n/(n-1)+2*sum(fit[[2]][,,k])/n
      }

    selmodel=order(AIC)[1]
    selindex=fit[[2]][,,selmodel]
    ord=order(selindex)
    selset=subset(ord,selindex[ord]!=0)

   return(list(selset,AIC,fit[[2]][,,selmodel]))
}

