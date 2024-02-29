###------------------------------------------------------------
#  The wrapper for Fortran 95 subroutines for analysis data
#   under the semiparametric odds ratio model
#   using the permutation likelihood approach
#   with some pre-fixed structure for parameters.
#
#  pmlkhfix is more general than pmlkh in two aspects.
#    a. allow for fixed zero parameters
#    b. compute and output log-likelihood
###---------------------------------------------------------------
#' The permutation likelihood approach with fixed structure
#'
#' This approach uses the maximum permutation likelihood approach
#' to estimating the parameters in the semiparametric odds ratio model
#' accommodating the pre-specified fixed structure in parameters.
#'
#' @param dat a data matrix of nxp dimension.
#' @param group a vector of positive integers of length ng
#'        such that sum(group)=p.
#' @param fixstruct a pxq matrix specifying location of fixed parameter (0),
#'                   or parameter to be estimated (1).
#' @param niter the maximum number of iterations.
#' @param eps convergence criterion: discrepancy between successive
#'            iterations for terminating the iteration.
#' @param nlag the averaged parameter values in the sequence to assess convergence, default nlag = 20.
#' @param nburnin the burnin steps before a sample is taken, default nburnin=5e4.
#' @param nsamp the sample sizes to be taken in the MCMC sampling of permutation, default nsamp = 1e4.
#' @param nintv the sampling interval (number of steps ignored) in the MCMC chain, default nintv = 2e2.
#' @param maxcyc the maximum step size (cycle of the permutation) in obtaining a candidate permutation
#'        (change from the previous one) in the Gibbs sampling, default maxcyc = 2.
#' @param steplen two-stage controls step length of the parameter update,steplen=1 means no control,
#'       steplen<1 means control.This parameter is set to make convergence more stable. the first
#'       component of steplen controls the step size in 1 to nstep iterations, the 2nd component
#'       of steplen controls the stepsize for the remaining iterations.
#' @param nstep controls the number of steps before a change of the step size occurs
#'
#' @details This method maximizes the permutation likelihood to obtain the parameter estimator
#'          and uses the inverse of the permutation information matrix
#'          to estimate the asymptotic variance of the estimator.
#'          The permutation likelihood is approximated using the Gibbs sampler
#'          of permutations with non-constant weights (probability).
#'
#' @return Estimate of the model parameters (column-wise vectorization) and
#'         the covariance matrix estimate corresponding to the vectorized parameters,
#'         and the log-likelihood.
#'
#' @references Chen, H.Y. (2022). Semiparametric Odds Ratio Model and its Application. CRC press.
#'
#' @examples \dontrun{
#' n=200; p=10
#' datmat=matrix(rnorm(n * p), ncol = p)
#' vargroup=c(2, 3, 1, 4)
#' structure <- matrix(rbinom(n=p*p,size=1,p=0.6), ncol = p)
#' pmlkhfix(dat=datmat, group=vargroup, fixstruct=structure)
#' }
#'
#' @export
pmlkhfix <- function(dat, group, fixstruct, niter = 50, eps = 1e-2, nlag = 20,
                     nburnin = 5e4, nsamp = 1e5, nintv = 5e1, maxcyc = 2,
                     steplen=c(1,1),nstep=15) {
  n <- dim(dat)[1]
  np <- dim(dat)[2]
  ng <- length(group)
  nq <- np * (np - 1)/2
  for(k in 1:ng){
    nq <- nq - group[k] * (group[k] - 1)/2
  } # nq=total number of parameters from blocks.

  theta <- rep(0, nq)
  estv <- matrix(0, nrow = nq, ncol = nq)

  theta2 <- array(0, c(niter, np, np))
  converge <- 0
  loglkh <- 0.0
  fit <- .Fortran("analysiswmgfix", as.double(dat), as.integer(n), as.integer(np),
                  as.integer(group), as.integer(ng), as.double(theta),
                  as.double(estv), as.integer(nq), as.double(eps),
                  as.integer(converge), as.double(loglkh), as.integer(fixstruct),
                  as.integer(niter), as.integer(nlag), as.integer(nburnin),
                  as.integer(nsamp), as.integer(nintv), as.integer(maxcyc),
                  as.double(theta2),as.double(steplen),as.integer(nstep))

  if(fit[[10]] == 0) {print("Convergence criterion is not met")}

  #print(fit[[6]])
  #print(matrix(fit[[7]],ncol=nq))
  draw(array(fit[[19]], c(niter, np, np)))

  return(list(fit[[6]], matrix(fit[[7]] , ncol = nq), fit[[11]], fit[[18]]))
}
