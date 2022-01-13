#' Pairwise pseudo-likelihood approach
#'
#' The \code{pwlkh} and \code{pwlkhfix} functions estimate the odds ratio parameters in the semiparametric odds ratio
#' model using the pairwise pseudo-likelihood approach proposed in Liang and Qin (2000).
#' The \code{pwpenlkh} function performs the node-wise selection of the network based on the semiparametric odds ratio
#' model. Currently, only the bilinear log-odds ratio function is allowed.
#'
#' @aliases pwlkh pwlkhfix pwpenlkh
#'
#' @param y the outcome data matrix. This must be specified.
#' @param x the covariate data matrix. This must be specified.
#' @param dat the whole data matrix.
#' @param lambda a vector of penalty parameters. A network is selected for each lambda value.
#' @param group a vector denoting the sizes of the variables in each group the data matrix divided.
#' @param fixstruct a matrix indicating which coefficients are set to 0.
#' @param niter the maximum number of iterations.
#' @param eps convergence criterion.
#' @param vect the method of vectoring parameter matrix. Either column-wise or row-wise.
#'
#' @usage pwlkh(y, x, niter = 50, eps = 1e-06, vect = "col")
#' @usage pwlkhfix(y, x, fixstruct, niter = 50, eps = 1e-06, vect = "col")
#' @usage pwpenlkh(dat, group, lambda, niter = 50, eps = 1e-06)
#'
#' @details
#' The \code{pwlkhfix} function is more general than the \code{pwlkh} function since the former not only allows for
#' the fixed zero coefficients, it also computes the logarithm of the pairwise pseudo-likelihood values.
#'
#' While the outcome data matrix and the variate data matrix has to be specified separately in the \code{pwlkh} and
#' the \code{pwlkhfix} functions, for the \code{pwpenlkh} function, one has to specify the whole data matrix including
#' the outcome and the covariates.
#'
#' For the \code{pwlkhfix} function, one has to specify a matrix for a fixed structure, that indicating which
#' coefficients are set to zero. For the \code{pwpenlkh} function, a \code{lambda} parameter that specifies the
#' penalties on the pairwise pseudo-likelihood, is required.
#'
#' @return
#' The \code{pwlkh} function returns
#'
#' 1. the odds ratio parameter estimates.
#'
#' 2. the estimated variance-covariance matrix of the odds ratio estimator.
#'
#' The \code{pwlkhfix} function returns
#'
#' 1. the odds ratio parameter estimates.
#'
#' 2. the estimated variance-covariance matrix of the odds ratio estimator.
#'
#' 3. the pairwise pseudo-likelihood value.
#'
#' The \code{pwpenlkh} function returns
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
#' n <- 200; p <- 3; q <- 4
#' covariates <- matrix(rnorm(n * q), ncol = q)
#' outcomes <- covariates %*% matrix(runif(q * p), ncol = p) + matrix(rnorm(n * p), ncol = p)
#'
#' # pwlkh
#' pwlkh(y = outcomes, x = covariates)
#'
#' # pwlkhfix
#' structure <- matrix(rbind(p * q), ncol = q)
#' pwlkhfix(y = outcomes, x = covariates, fixstruct = structure)
#'
#' # pwpenlkh
#' n <- 200; p <- 10
#' datmat <- matrix(rnorm(n * p), ncol = p)
#' vargroup <- c(2, 3, 1, 4)
#' penaltyparam <- c(10, 50, 500)
#' pwpenlkh(dat = datmat, group = vargroup, lambda = penaltyparam)
#'
#' @references Liang, K. Y. and Qin, J. (2000). Regression analysis under non-standard situations: a
#'   pseudolikelihood approach. Journal of the Royal Statistical Society, Ser. B, 62:773-786.
#' @references Chen, H. Y. (2007). A semiparametric odds ratio model for measuring association.
#'    Biometrics, 63: 413-421.
#' @references ...
#'
# Simple approach
#' @export
pwlkh <- function(y, x, niter = 50, eps = 1e-6, vect = "col") {

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
                    as.integer(niter), as.double(eps), as.integer(converge))
  } else {
    fit <- .Fortran("pwmlerow", as.double(y), as.double(x), as.integer(n), as.integer(np),
                    as.integer(nq), as.double(theta), as.double(estv),
                    as.integer(niter), as.double(eps), as.integer(converge))
  }

  if(fit[[10]] == 0) {print("Convergence criterion is not met")}

  return(list(fit[[6]], matrix(fit[[7]], ncol = np * nq)))
}
# Pairwise pseudo-likelihood approach with fixed structure
#' @export
pwlkhfix <- function(y, x, fixstruct, niter=50, eps=1e-6, vect="col") {

  if(is.matrix(y) == TRUE) {
    n <- dim(y)[1]
    np <- dim(y)[2]
  } else{
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
  loglkh <- 0
  if(vect == "col"){
    fit <- .Fortran("pwmlecolfix", as.double(y), as.double(x), as.integer(n), as.integer(np),
                    as.integer(nq), as.integer(fixstruct),
                    as.double(theta), as.double(estv), as.double(loglkh),
                    as.integer(niter), as.double(eps), as.integer(converge))
  } else{
    fit <- .Fortran("pwmlerowfix", as.double(y), as.double(x), as.integer(n), as.integer(np),
                    as.integer(nq), as.integer(fixstruct),
                    as.double(theta), as.double(estv), as.double(loglkh),
                    as.integer(niter), as.double(eps), as.integer(converge))
  }
  if(fit[[12]] == 0) {print("Convergence criterion is not met")}

  return(list(fit[[7]], matrix(fit[[8]], ncol = np * nq), fit[[9]]))
}
# Penalized pseudo-likelihood approach for node-wise selection
#' @export
pwpenlkh <- function(dat, group, lambda, niter = 50, eps = 1e-6) {

  n <- dim(dat)[1]
  node <- dim(dat)[2]
  ng <- length(group)
  nlam <- length(lambda)

  selnet <- array(0, c(node,node,nlam))

  fit <- .Fortran("networkselectbypw", as.double(dat), as.integer(n), as.integer(node),
                  as.integer(group), as.integer(ng), as.double(lambda), as.integer(nlam),
                  as.integer(niter), as.double(eps), as.integer(selnet))

  #create network graph#
  network <- array(fit[[10]], c(node,node,nlam))
  networka <- network

  for(j in 1:nlam){
    ncurrent <- 0
    for(k in 1:ng) {# blockout the within group pairs
      network[ncurrent + c(1:group[k]), ncurrent + c(1:group[k]), j] = matrix(NA, nrow = group[k], ncol = group[k])
      ncurrent = ncurrent + group[k]
    }
  }
  annotation <- c("0 = Not connected", "1 = weak connection (detected once)", "2 = strong connection (detected twice)")

  return(list(group, network, annotation, networka))
}

