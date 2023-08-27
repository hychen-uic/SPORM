###------------------------------------------------------------
#  The wrapper for Fortran 95 subroutines for analysis data
#      under the semiparametric odds ratio model
#      using the semiparametric likelihood approach
###---------------------------------------------------------------
#' The simple semiparametric likelihood approach
#'
#' This approach uses the maximum semiparametric likelihood approach to estimating
#' the parameters in the semiparametric odds ratio model.
#'
#' @param y outcomes: can be a vector of length n or a matrix of nxq dimension.
#' @param x covariates: a matrix of nxp dimension.
#' @param niter maximum number of iterations in finding the estimator.
#' @param eps convergence criterion: discrepancy between successive
#'            iterations for terminating the iteration.
#' @param vect methods of vectorization of parameter matrix.
#'        Either column-wise (vect='col') or row-wise (vect='row').
#'
#' @details This method maximizes the semiparametric likelihood to obtain the parameter estimator
#'          and uses the inverse of the profile information matrix to estimate
#'          the asymptotic variance of the estimator.
#'
#' @return Estimate of the model parameters (vectorized) and
#'         the covariance matrix estimate corresponding to the vectorized parameters.
#'
#' @references Chen, H. Y. (2007). A semiparametric odds ratio model for measuring association,
#'             *Biometrics*,**63**,413--421.
#' @references Chen, H.Y. (2022). Semiparametric Odds Ratio Model and its Application. CRC press.
#' @references Chen, H. Y., and Rader, D. E., and Li, M. (2015). Likelihood inference on
#'              semiparametric odds ratio Model. *Journal of the American Statistical Association*,
#'              **110**,1125-1135.
#'
#' @examples \dontrun{Example 1:
#' n=100; p=2
#' x=matrix(2*runif(n*p)-1,ncol=p)
#' y=x[,1]-x[,2]+rnorm(n)
#' splkh(y,x,niter=50,eps=1e-6,vect='col')
#'
#' Example 2:
#' n=200;p=2;q=2
#' x=matrix(2*runif(n*p)-1,ncol=p)
#' y=array(0,c(n,q))
#' y[,1]=x[,1]-x[,2]+rnorm(n)
#' y[,2]=x[,1]+x[,2]+rnorm(n)
#' splkh(y,x,niter=40,eps=1e-6,vect="row")}
#'
#' @export
splkh <- function(y, x, niter = 50, eps = 1e-6, vect = "col") {

  if (is.matrix(y) == TRUE) {
    n <- dim(y)[1]
    np <- dim(y)[2]
  } else{
     n <- length(y)
     np <- 1
  }
  if (is.matrix(x) == TRUE) {
    nq <- dim(x)[2]
  } else {
     nq <- 1
  }

  theta <- rep(0, np * nq)
  estv <- matrix(0, nrow = np * nq, ncol = np * nq)

  converge <- 0
  if (vect == "col") {
    fit <- .Fortran("ormlecol", as.double(y), as.double(x), as.integer(n), as.integer(np),
                 as.integer(nq), as.double(theta), as.double(estv),
                 as.integer(niter), as.double(eps), as.integer(converge))
  } else {
    fit <- .Fortran("ormlerow", as.double(y), as.double(x), as.integer(n), as.integer(np),
                 as.integer(nq), as.double(theta), as.double(estv),
                 as.integer(niter), as.double(eps), as.integer(converge))
  }

  if(fit[[10]] == 0) {print("Convergence criterion is not met")}

  return(list(fit[[6]], matrix(fit[[7]], ncol = np * nq)))
}
