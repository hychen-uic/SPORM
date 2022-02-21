###------------------------------------------------------------
#  The wrapper for Fortran 95 subroutines for network selection
#    under the semiparametric odds ratio model
#    using the penalized semiparametric likelihood approach.
#    The network selection is the node-wise approach.
###---------------------------------------------------------------
#' The penalized semiparametric likelihood approach for node-wise network selection
#'
#' This approach selects network structure
#'    under the semiparametric odds ratio model
#'    using the penalized semiparametric likelihood approach.
#'    The network selection is the node-wise approach.
#'
#' @param dat a data matrix of nxnode dimension
#' @param group a vector of positive integers of length ng
#'        such that sum(group)=node
#' @param lambda a vector of penalty values for network selection,
#'               each penalty value determine a network
#' @param niter maximum number of iterations in finding the estimator
#' @param eps convergence criterion: discrepancy between successive
#'            iterations for terminal the iteration
#'
#' @details This method maximizes the penalized semiparametric likelihood
#'          with penalty being the l^1-norm of the parameters to select
#'          the network structure.
#' @return group denoting the input node clusters, a set of networks identified corresponding
#'           to the penalty values (with within group connections coded by NA),
#'           annotation of network connection (0=not connected etc.)
#'           the set of networks identified with within group connections coded by 0.
#'
#' @references Chen, H.Y. (2022). Semiparametric Odds Ratio Model and its Application. CRC press.
#' @references Chen, H. Y. and Chen, J. (2020). Network selection through semiparametric odds ratio model. Manuscript.
#'
#' @examples \dontrun{
#' Use the internal data file name dat, a 400x9 data matrix.
#' group=c(2,3,4)
#' lambda=c(10,50,200,500,800)
#' sppenlkh(dat=dat,group=group,lambda=lambda)
#' }
#'
#' @export
sppenlkh <- function(dat, group, lambda, niter = 50, eps = 1e-6) {

  n <- dim(dat)[1]
  node <- dim(dat)[2]
  ng <- length(group)
  nlam <- length(lambda)

  selnet<- array(0, c(node, node, nlam))

  fit <- .Fortran("networkselectbysp", as.double(dat), as.integer(n), as.integer(node),
                  as.integer(group), as.integer(ng), as.double(lambda), as.integer(nlam),
                  as.integer(niter), as.double(eps), as.integer(selnet))


  #create network graph#
  network <- array(fit[[10]], c(node,node,nlam))
  networka <- network

  for(j in 1:nlam){
    ncurrent <- 0
    for(k in 1:ng){ #blockout the within group pairs
      network[ncurrent + c(1:group[k]), ncurrent + c(1:group[k]),j] = matrix(NA, nrow = group[k],ncol = group[k])
      ncurrent <- ncurrent + group[k]
    }
  }
  annotation <- c("0 = Not connected", "1 = weak connection (detected once)", "2 = strong connection (detected twice)")

  return(list(group, network, annotation, networka))
}

