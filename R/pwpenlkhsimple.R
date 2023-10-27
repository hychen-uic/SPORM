###---------------------------------------------------------------
#  The wrapper for Fortran 95 subroutines for model selection
#    under the semiparametric odds ratio model
#    using the penalized pairwise pseudo-likelihood approach.
###---------------------------------------------------------------
#' The penalized pairwise pseudo-likelihood approach for model selection
#'
#' This approach selects a model by parameter selection
#'    under the semiparametric odds ratio model
#'    using the penalized pairwise pseudo-likelihood approach.
#'
#' @param y a vector of outcomes
#' @param x a vector of covariates
#' @param lambda a vector of penalty values for model selection,
#'               each penalty value determine a model.
#' @param niter maximum number of iterations in finding the estimator.
#' @param eps convergence criterion: discrepancy between successive
#'            iterations for terminal the iteration
#'
#' @details This method maximizes the penalized pairwise pseudo-likelihood
#'          with penalty being the l^1-norm of the parameters to select
#'          a sparse model.
#' @return a set of selected parameters to determine the model.
#'
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
#' pwpenlkhsimple(y,x,lambda=lambda)
#' }
#'
#'@export
pwpenlkhsimple <- function(y,x,lambda,niter=50,eps=1e-6) {
  if(!is.loaded("pairwisePENlkh.so")){
    dyn.load("src/makevars/pairwisePENlkh.so")
  }else{
    dyn.unload("src/makevars/pairwisePENlkh.so")
    dyn.load("src/makevars/pairwisePENlkh.so")
  }
  n=dim(y)[1]
  p=dim(y)[2]
  q=dim(x)[2]
  nlam=length(lambda)
  theta=array(0,c(p,q))
  selparm=array(0,c(p,q,nlam))

  for(k in 1:nlam){
    fit=.Fortran("penpwlkh",as.double(y),as.double(x),
                  as.integer(n),as.integer(p),as.integer(q),
                  as.double(theta),as.integer(selparm[,,k]),
                  as.integer(niter),as.double(eps),as.double(lambda[k]))
    }
  dyn.unload("src/makevars/pairwisePENlkh.so")

  return(list(lambda,selparm))
}
