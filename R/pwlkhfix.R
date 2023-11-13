###------------------------------------------------------------
#  The wrapper for Fortran 95 subroutines for analysis data
#   under the semiparametric odds ratio model
#   using the pairwise pseudo-likelihood approach
#   with some pre-fixed structure for parameters.
#
#  pwlkhfix is more general than splkh in two aspects.
#    a. allow for fixed zero parameters
#    b. compute and output log-likelihood
###---------------------------------------------------------------
#' The pairwise pseudo-likelihood approach with fixed structure
#'
#' This approach uses the maximum pairwise pseudo-likelihood approach
#' to estimating the parameters in the semiparametric odds ratio model
#' accommodating the pre-specified fixed structure in parameters.
#'
#' @param y outcomes: a matrix of nxq dimension.
#' @param x covariates: a matrix of nxp dimension.
#' @param theta a pxq matrix specifying fixed parameter values. non-fixed values do not matter
#'        the default value theta=0 is for convenience and is automatically redefined appropriately
#'        For non-default values, theta needs to be conformed with the actual theta in dimension.
#' @param fixstruct a pxq matrix specifying location of parameter to be fixed (0)
#'                   or parameter to be estimated (1).
#'                   NOTE: Fixed parametervalues can be nonzero and determined by theta.
#' @param niter maximum number of iterations in finding the estimator.
#' @param eps convergence criterion: discrepancy between successive
#'            iterations for terminating the iteration.
#' @param vect methods of vectorization of parameter matrix.
#'        Either column-wise (vect='col') or row-wise (vect='row').
#'
#' @details This method maximizes the pairwise likelihood to obtain the parameter estimator
#'          and uses U-statistic theory to estimate the asymptotic variance of the estimator.
#'
#' @return Estimate of the model parameters (vectorized),
#'         the covariance matrix estimate corresponding to the vectorized parameters,
#'         and the log-likelihood.
#'
#' @references Chan K. C. G.(2013). Nuisance parameter elimination for proportional likelihood
#'             ratio models with nonignorable missingness and random truncation. Biometrika, 100,
#'             269-276.
#' @references Chen, H.Y. (2022). Semiparametric Odds Ratio Model and its Application. CRC press.
#'
#' @references Liang, K.Y. and Qin, J. (2000). Regression analysis under non-standard situations:
#'             a pairwise pseudolikelihood approach, *Journal of the Royal Statistical Society, Ser. B*,
#'             **62**, 773-786.
#'
#' @examples \dontrun{
#' # Use the internal data file name dat, a 400x9 data matrix.
#' outcomes=dat[,1:2]
#' covariates=dat[,4:7]
#' structure <- matrix(rbinom(n=2*4,size=1,p=0.6), ncol = 4)
#' theta=matrix(0,nrow=2,ncol=4)
#' pwlkhfix(y = outcomes, x = covariates, fixstruct = structure)
#' pwlkhfix(y,x,fixstruct=structure,theta=c(rep(0.1,4),rep(0,4)))
#' }
#'
#' @export
#'
pwlkhfix <- function(y, x, fixstruct, theta=0, niter=50, eps=1e-6, vect="col") {

  if(is.matrix(y) == TRUE) {
    n <- dim(y)[1]
    np <- dim(y)[2]
  } else{
    n <- length(y)
    np <- 1
  }
  if (is.matrix(x) == TRUE) {
    nq <- dim(x)[2]
  }else{
    nq <- 1
  }
  if(sum(abs(theta))<1e-8){ # allow default value
    theta <- rep(0, np * nq)
    }
  estv <- matrix(0, nrow = np * nq, ncol = np * nq)

  converge <- 0
  loglkh <- 0
  if(vect == "col"){
    fit <- .Fortran("pwmlecolfix",as.double(y), as.double(x), as.integer(n), as.integer(np),
                    as.integer(nq), as.integer(fixstruct),
                    as.double(theta), as.double(estv), as.double(loglkh),
                    as.integer(niter), as.double(eps), as.integer(converge))
  }else{
    fit <- .Fortran("pwmlerowfix", as.double(y), as.double(x), as.integer(n), as.integer(np),
                    as.integer(nq), as.integer(fixstruct),
                    as.double(theta), as.double(estv), as.double(loglkh),
                    as.integer(niter), as.double(eps), as.integer(converge))
  }
  if(fit[[12]] == 0) {print("Convergence criterion is not met")}

  return(list(fit[[7]], matrix(fit[[8]], ncol = np * nq), fit[[9]]))
}
