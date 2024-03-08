###------------------------------------------------------------
#  The wrapper for Fortran 95 subroutines for analysis data
#      under the semiparametric odds ratio model
#      using the pairwise pseudo-likelihood approach
###---------------------------------------------------------------
#' The simple pairwise pseudo-likelihood approach
#'
#' This approach uses the maximum pairwise pseudo-likelihood approach to estimating
#' the parameters in the semiparametric odds ratio model.
#'
#' @param y outcomes: can be a vector of length n or a matrix of nxq dimension.
#' @param x covariates: a matrix of nxp dimension.
#' @param niter maximum number of iterations in finding the estimator.
#' @param eps convergence criterion: discrepancy between successive
#'            iterations for terminating the iteration.
#' @param vect methods of vectorization of parameter matrix.
#'        Either column-wise (vect='col') or row-wise (vect='row').
#' @param steplen two-stage controls step length of the parameter update,steplen=1 means no control,
#'       steplen<1 means control.This parameter is set to make convergence more stable. the first
#'       component of steplen controls the step size in 1 to nstep iterations, the 2nd component
#'       of steplen controls the stepsize for the remaining iterations.
#' @param nstep controls the number of steps before a change of the step size occurs
#'
#' @details This method maximizes the pairwise likelihood to obtain the parameter estimator
#'          and uses U-statistic theory to estimate the asymptotic variance of the estimator.
#'
#' @return Estimate of the model parameters (vectorized) and
#'         the covariance matrix estimate corresponding to the vectorized parameters.
#'
#' @references Chan K. C. G. (2013). Nuisance parameter elimination for proportional likelihood
#'             ratio models with nonignorable missingness and random truncation. *Biometrika*, **100**,
#'             269-276.
#'
#' @references Chen, H.Y. (2022). Semiparametric Odds Ratio Model and its Application. CRC press.
#'
#' @references Liang, K.Y. and Qin, J. (2000). Regression analysis under non-standard situations:
#'             a pairwise pseudolikelihood approach. *Journal of the Royal Statistical Society, Ser. B*,
#'             **62**, 773-786.
#'
#' @examples \dontrun{
#' Example 1:
#' n=100; p=2
#' x=matrix(2*runif(n*p)-1,ncol=p)
#' y=x[,1]-x[,2]+rnorm(n)
#' pwlkh(y,x,niter=50,eps=1e-6,vect='col')
#'
#' Example 2:
#' n=200;p=2;q=2
#' x=matrix(2*runif(n*p)-1,ncol=p)
#' y=array(0,c(n,q))
#' y[,1]=x[,1]-x[,2]+rnorm(n)
#' y[,2]=x[,1]+x[,2]+rnorm(n)
#' pwlkh(y,x,niter=40,eps=1e-6,vect="row")
#' }
#'
#' @export
pwlkh <- function(y, x, niter = 50, eps = 1e-6, vect = "col",steplen=c(0.5,1),nstep=10) {

  if (is.matrix(y) == TRUE) {
    n <- dim(y)[1]
    np <- dim(y)[2]
  } else {
    n <- length(y)
    np <- 1
  }
  if (is.matrix(x) == TRUE) {
    nq <- dim(x)[2]
  } else{
    nq <- 1
  }

  theta <- rep(0, np * nq)
  estv <- matrix(0, nrow = np * nq, ncol = np * nq)

  converge <- 0

  if (vect == "col"){
    fit <- .Fortran("pwmlecol", as.double(y), as.double(x), as.integer(n), as.integer(np),
                    as.integer(nq), as.double(theta), as.double(estv),
                    as.integer(niter), as.double(eps), as.integer(converge),as.double(steplen), as.integer(nstep))
  } else {
    fit <- .Fortran("pwmlerow", as.double(y), as.double(x), as.integer(n), as.integer(np),
                    as.integer(nq), as.double(theta), as.double(estv),
                    as.integer(niter), as.double(eps), as.integer(converge),as.double(steplen), as.integer(nstep))
  }

  if(fit[[10]] == 0) {print("Convergence criterion is not met")}

  return(list(fit[[6]], matrix(fit[[7]], ncol = np * nq)))
}
