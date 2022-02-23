###------------------------------------------------------------
#  The wrapper for Fortran 95 subroutines for analysis data
#   under the semiparametric odds ratio model
#   using the semiparametric likelihood approach
#   with some pre-fixed structure for parameters.
#
#  splkhfix is more general than splkh in two aspects.
#    a. allow for fixed zero parameters
#    b. compute and output log-likelihood
###---------------------------------------------------------------
#' The semiparametric likelihood approach with fixed structure
#'
#' This approach uses the maximum semiparametric likelihood approach
#' to estimating the parameters in the semiparametric odds ratio model
#' accommodating the pre-specified fixed structure in parameters.
#'
#' @param y outcomes: a matrix of nxq dimension.
#' @param x covariates: a matrix of nxp dimension.
#' @param fixstruct a pxq matrix specifying location of fixed parameter (0),
#'                   or parameter to be estimated (1).
#' @param niter the maximum number of iterations in finding the estimator.
#' @param eps convergence criterion: discrepancy between successive
#'            iterations for terminating the iteration.
#' @param vect methods of vectorization of parameter matrix.
#'        Either column-wise (vect='col') or row-wise (vect='row').
#'
#' @details This method maximizes the semiparametric likelihood to obtain the parameter estimator
#'          and uses the inverse of the profile information matrix
#'          to estimate the asymptotic variance of the estimator.
#'
#' @return Estimate of the model parameters (vectorized),
#'         the covariance matrix estimate corresponding to the vectorized parameters,
#'         and the log-likelihood.
#'
#' @references Chen, H.Y. (2022). Semiparametric Odds Ratio Model and its Application. CRC press.
#' @references Chen, H. Y., and Rader, D. E., and Li, M. (2015). Likelihood inference on
#'              semiparametric odds ratio Model. *Journal of the American Statistical Association*,
#'              **110**,1125-1135.
#'
#' @examples \dontrun{
#' # Use the internal data file name dat, a 400x9 data matrix.
#' outcomes=dat[,1:2]
#' covariates=dat[,4:7]
#' structure <- matrix(rbinom(n=2*4,size=1,p=0.6), ncol = 4)
#' splkhfix(y = outcomes, x = covariates, fixstruct = structure)
#' }
#'
#' @export
splkhfix <- function(y, x, fixstruct, niter=50, eps=1e-6, vect="col") {

  if(is.matrix(y) == TRUE) {
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
  loglkh <- 0
  if(vect == "col") {
    fit <- .Fortran("ormlecolfix", as.double(y), as.double(x), as.integer(n), as.integer(np),
                    as.integer(nq), as.integer(fixstruct),
                    as.double(theta), as.double(estv), as.double(loglkh),
                    as.integer(niter), as.double(eps), as.integer(converge))
  }
  else{
    fit <- .Fortran("orMLErowfix", as.double(y), as.double(x), as.integer(n), as.integer(np),
                    as.integer(nq), as.integer(fixstruct),
                    as.double(theta), as.double(estv), as.double(loglkh),
                    as.integer(niter), as.double(eps), as.integer(converge))
  }

  if(fit[[12]] == 0) {print("Convergence criterion is not met")}

  return(list(fit[[7]], matrix(fit[[8]], ncol = np * nq), fit[[9]]))
}
