###------------------------------------------------------------
#  The wrapper for Fortran 95 subroutines for analysis data
#      under the semiparametric odds ratio model
#      using the permutation likelihood approach
#      through the Monte Carlo sampling approximation
###---------------------------------------------------------------
#' The simple permutation likelihood approach
#'
#' This approach uses the maximum permutation likelihood approach to estimating
#' the parameters in the semiparametric odds ratio model.
#'
#' @param dat a data matrix of nxp dimension.
#' @param group a vector of positive integers of length ng
#'        such that sum(group)=p.
#' @param niter the maximum number of iterations.
#' @param eps convergence criterion: discrepancy between successive
#'            iterations for terminating the iteration.
#' @param nlag the averaged parameter values in the sequence to assess convergence, default nlag=20.
#' @param plot Default is TRUE, to produce the convergence plot. If FALSE, suppress the plot.
#' @param nburnin the burnin steps before a sample is taken, default nburnin=5e4.
#' @param nsamp the sample sizes to be taken in the MCMC sampling of permutation, default nsamp=1e4.
#' @param nintv the sampling interval (number of steps ignored) in the MCMC chain, default nintv=2e2.
#' @param maxcyc the maximum step size (cycle of the permutation) in obtaining a candidate permutation
#'        (change from the previous one) in the Gibbs sampling, default maxcyc=2.
#'
#' @details This method maximizes the permutation likelihood to obtain the parameter estimator
#'          and uses the inverse of the permutation likelihood information matrix to estimate
#'          the asymptotic variance of the estimator.
#'          The permutation likelihood is approximated using the Gibbs sampler
#'          of permutations with non-constant weights (probability).
#'
#' @return Estimate of the model parameters (column-wise vectorization) and
#'         the covariance matrix estimate corresponding to the vectorized parameters.
#'
#' @references Chen, H.Y. (2022). Semiparametric Odds Ratio Model and its Application. CRC press.
#'
#' @examples \dontrun{
#' n=200; p=10
#' datmat=matrix(rnorm(n * p), ncol = p)
#' vargroup=c(2, 3, 1, 4)
#' pmlkh(dat=datmat, group=vargroup)
#' }
#'
#' @export
pmlkh <- function(dat, group, niter = 50, eps = 1e-2,nlag = 20, plot=TRUE,
                  nburnin = 5e4, nsamp = 1e5, nintv = 5e1, maxcyc = 2) {
  n <- dim(dat)[1]
  np <- dim(dat)[2]
  ng <- length(group)
  nq <- np * (np-1)/2
  for(k in 1:ng){
    nq <- nq - group[k] * (group[k] - 1)/2
  } # nq=total number of parameters from blocks.

  theta <- rep(0, nq)
  estv <- matrix(0, nrow = nq, ncol = nq)

  theta2 <- array(0, c(niter, np, np))
  converge <- 0

  fit <- .Fortran("analysiswmg", as.double(dat), as.integer(n), as.integer(np),
               as.integer(group), as.integer(ng), as.double(theta),
               as.double(estv), as.integer(nq), as.double(eps),
               as.integer(converge), as.integer(niter), as.integer(nlag),
               as.integer(nburnin), as.integer(nsamp), as.integer(nintv),
               as.integer(maxcyc), as.double(theta2))

  if(fit[[10]] == 0){print("Convergence criterion is not met")}
  #print(fit[[6]])
  #print(matrix(fit[[7]],ncol=nq))
  #print(array(fit[[17]],c(niter,np,np)))
  if(plot==TRUE){
    draw(array(fit[[17]], c(niter, np, np)))
    }
  return(list(fit[[6]], matrix(fit[[7]], ncol = nq), fit[[16]]))
}

