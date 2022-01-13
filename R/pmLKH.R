#' Permutation likelihood approach
#'
#' The \code{pmlkh} and \code{pmlkhfix} functions maximize the Monte Carlo approximated permutation likelihood
#' to obtain the parameter estimates in the semiparametric odds ratio model.  The \code{pmpenlkh} function
#' uses the penalized approximated permutation likelihood approach to select network based on the semiparametric
#' odds ratio model. Currently, only the bilinear log-odds ratio function is allowed.
#'
#' @aliases pmlkh pmlkhfix pmpenlkh
#'
#' @param dat the whole data matrix.
#' @param group a vector denoting the sizes of the variables in each group the data matrix divided.
#' @param lambda a vector of penalty parameters. A network is selected for each lambda value.
#' @param fixstruct a matrix indicating which coefficients are set to 0.
#' @param niter the maximum number of iterations.
#' @param eps convergence criterion.
#' @param nlag the number of steps averaged in obtaining the stablized results.
#' @param nburnin samples to be discarded in the metropolis algorithm for sampling permutations.
#' @param nsamp size of the Monte Carlo sample taken for mean and variance matrix approximation.
#' @param nintv sampling interval in the Monte Carlo Markov chain (Metropolis).
#' @param maxcyc the maximum cycle length in the candidate permutation update.
#'
#' @usage pmlkh(dat, group, niter = 50, eps = 0.01, nlag = 20, nburnin = 50000,
#'       nsamp = 1e+05, nintv = 50, maxcyc = 2)
#' @usage pmlkhfix(dat, group, fixstruct, niter = 50, eps = 0.01, nlag = 20, nburnin = 50000,
#'          nsamp = 1e+05, nintv = 50, maxcyc = 2)
#' @usage pmpenlkh(dat, group, lambda, nburnin = 50000, nsamp = 10000,
#'          nintv = 200, niter = 50, eps = 0.01, maxcyc = 2)
#'
#' @details
#' The \code{pmlkhfix} function is more general than the \code{pmlkh} function since the former not only allows for
#' the fixed zero coefficients, it also computes the logarithm of the maximum permutation likelihood values.
#'
#' For all of the three functions, one has to specify the whole data matrix including the outcome and covariates,
#' and the sizes of the variables in each group the data matrix divides.
#'
#' For the \code{pmlkhfix} function, one has to specify a matrix for a fixed structure, that indicating which
#' coefficients are set to zero. For the \code{pmpenlkh} function, a \code{lambda} parameter that specifies the
#' penalies on the permutation likelihood, is required.
#'
#' While the \code{pmlkh} and \code{pmlkhfix} are mainly used for estimating the odds ratio parameters and the
#' associated variance-covariance matrices, the \code{pmpenlkh} function performs the selection of network
#' connections.
#'
#' @return
#' The \code{pmlkh} function returns
#'
#' 1. the estimated odds ratio parameters.
#'
#' 2. the variance-covariance matrix of the odds ratio parameter estimator.
#'
#' The \code{pmlkhfix} function returns
#'
#' 1. the estimated odds ratio parameters.
#'
#' 2. the variance-covariance matrix of the odds ratio parameter estimator.
#'
#' 3. logarithm of the maximum permutation likelihood.
#'
#' The \code{pmpenlkh} function returns
#'
#' 1. The group vector specified.
#'
#' 2. A set of networks selected by the penalized pseudo-likelihood with within group labeled by NA.
#'
#' 3. Annotation of the entries of the network structure matrix.
#'
#' 4. A set of networks selected by the penalized pseudo-likelihood with within group labeled by 0.
#'
#' @examples
#'
#' \dontrun{
#' # pmklh example 1
#' n <- 100
#' p <- 4
#' group <- c(2, 2)
#' dat <- array(0, c(n, p))
#' dat[, (1 + group[1]) : p] <- matrix((runif(n * group[2]) - 0.5) * 2, ncol = group[2])
#' dat[, 1] <- rpois(n, exp(dat[,3] - dat[,4])) #Poisson error
#' dat[, 2] <- rpois(n, exp(dat[,3] + dat[,4]))
#'
#' pmlkh(dat, group)
#'
#' # pmlkh example 2
#' n <-  200
#' p <- 9
#' group <- c(2, 3, 4)
#' dat <- matrix(rnorm(n * p), ncol = p)
#' for(k in 1:p - 2){
#'   dat[, k] <- dat[, k + 2] + rnorm(n)
#' }
#'
#' pmlkh(dat, group, 60, 30, 40000, 5000, 300)
#'
#' # pmlkhfix example 1
#' n <- 100
#' p <- 4
#' vargroup <- c(2, 2)
#' datmat <- array(0, c(n, p))
#' datmat[, (1+group[1]) : p] <- matrix((runif(n * group[2]) - 0.5) * 2, ncol = group[2])
#' datmat[,1] <- rpois(n,exp(dat[,3] - dat[,4])) #Poisson error
#' datmat[,2] <- rpois(n,exp(dat[,3] + dat[,4]))
#' structure <- matrix(rbin(p*p),ncol=p)
#'
#' pmlkhfix(dat = datmat, group = vargroup, fixstruct = structure)
#'
#' # pmlkhfix example 2
#' n <- 200
#' p <- 9
#' vargroup <- c(2, 3, 4)
#' datmat <- matrix(rnorm(n * p), ncol = p)
#' structure <- matrix(bin(p*p),ncol=p)
#' for(k in 1:p - 2){
#'   datmat[,k] <- datmat[,k+2] + rnorm(n)
#' }
#'
#' pmlkhfix(dat = datmat, group = vargroup, fixstruct = structure,
#'          niter = 60, nlag = 30, nburn = 40000, nsamp = 5000, nintv = 300)
#'
#' # pmpenlkh example
#' n <- 200; p <- 10
#' datmat <- matrix(rnorm(n * p), ncol = p)
#' vargroup <- c(2, 3, 1, 4)
#' penaltyparam <- c(10, 50, 500)
#'
#' pmpenlkh(dat = datmat, group = vargroup, lambda = penaltyparam)
#'
#' }
#'
#' @references ...
#' @references ...
#'
# Simple approach
#' @export
pmlkh <- function(dat, group, niter = 50, eps = 1e-2, nlag = 20,
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
  draw(array(fit[[17]], c(niter, np, np)))
  return(list(fit[[6]], matrix(fit[[7]], ncol = nq), fit[[16]]))
}
# Permutation likelihood approach with fixed structure
#' @export
pmlkhfix <- function(dat, group, fixstruct, niter = 50, eps = 1e-2, nlag = 20,
                     nburnin = 5e4, nsamp = 1e5, nintv = 5e1, maxcyc = 2) {
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
               as.double(theta2))

  if(fit[[10]] == 0) {print("Convergence criterion is not met")}

  #print(fit[[6]])
  #print(matrix(fit[[7]],ncol=nq))
  draw(array(fit[[19]], c(niter, np, np)))

  return(list(fit[[6]], matrix(fit[[7]] , ncol = nq), fit[[11]], fit[[18]]))
}
# Penalized permutation likelihood approach for selecting networks
#' @export
pmpenlkh <- function(dat, group, lambda, nburnin = 5e4,
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
               as.integer(nburnin), as.integer(nsamp), as.integer(nintv),
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
