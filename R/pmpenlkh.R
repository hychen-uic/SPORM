###------------------------------------------------------------
#  The wrapper for Fortran 95 subroutines for network selection
#    under the semiparametric odds ratio model
#    using the penalized permutation likelihood approach.
#    The network selection is the node-wise approach.
###---------------------------------------------------------------
#' The penalized permutation likelihood approach for node-wise network selection
#'
#' This approach selects network structure
#'    under the semiparametric odds ratio model
#'    using the penalized permutation likelihood approach.
#'    The network selection is the joint selection.
#'
#' @param dat a data matrix of nxnode dimension.
#' @param group a vector of positive integers of length ng
#'        such that sum(group)=node.
#' @param lambda a vector of penalty values for network selection,
#'               each penalty value determines a network.
#' @param burnin samples to be discarded in the Metropolis algorithm for sampling permutations.
#' @param nsamp size of the Monte Carlo sample taken for mean and variance matrix approximation.
#' @param nintv sampling interval in the Monte Carlo Markov chain (Metropolis).
#' @param niter maximum number of iterations in finding the estimator.
#' @param eps convergence criterion: discrepancy between successive
#'            iterations for terminating the iteration
#' @param maxcyc the maximum cycle length in the candidate permutation update.
#'
#' @details This method maximizes the penalized permutation likelihood
#'          with penalty being the l^1-norm of the parameters to select
#'          the network structure. The selection is joint selection,
#'          NOT a node-wise selection, and is therefore more efficient statistically.
#'
#' @return Group denoting the input node clusters, a set of networks identified corresponding
#'           to the penalty values (with within group connections coded by NA),
#'           annotation of network connection (0=not connected etc.)
#'           the set of networks identified with within group connections coded by 0.
#'
#' @references Chen, H.Y. (2022). Semiparametric Odds Ratio Model and its Application. CRC press.
#' @references Chen, H. Y. and Chen, J. (2020). Network selection through semiparametric odds ratio model. Manuscript.
#'
#' @examples \dontrun{
#' n=200; p=10
#' datmat=matrix(rnorm(n * p), ncol = p)
#' vargroup=c(2, 3, 1, 4)
#' penaltyparam=c(10, 50, 500)
#' pmpenlkh(dat=datmat, group=vargroup, lambda=penaltyparam)
#'}
#'
#' @export
pmpenlkh <- function(dat, group, lambda, burnin = 5e4,
                     nsamp = 1e4,nintv = 2e2, niter = 50, eps = 1e-2, maxcyc = 2) {

  n <- dim(dat)[1]
  node <- dim(dat)[2]
  ng <- length(group)
  nlam <- length(lambda)

  selnet <- array(0, c(node, node, nlam))

  fit <- .Fortran("networkselectbypm", as.double(dat), as.integer(n),
                  as.integer(node), as.integer(group), as.integer(ng),
                  as.double(lambda), as.integer(nlam),
                  as.integer(niter), as.double(eps), as.integer(selnet),
                  as.integer(burnin), as.integer(nsamp), as.integer(nintv),
                  as.integer(maxcyc))

  #create network graph#
  network <- array(fit[[10]], c(node, node, nlam))
  networka <- network

  for (j in 1:nlam) {
    ncurrent = 0
    for (k in 1:ng) {#blockout the within group paairs
      network[ncurrent + c(1:group[k]), ncurrent + c(1:group[k]), j] = matrix(NA, nrow = group[k], ncol = group[k])
      ncurrent = ncurrent + group[k]
    }
  }
  annotation <- c("0 = Not connected", "1 = weak connection (detected once)", "2 = strong connection (detected twice)")

  return(list(group, network, annotation, networka))
}
