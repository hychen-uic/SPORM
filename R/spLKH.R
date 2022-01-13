#' Semiparametric likelihood approach
#'
#' The \code{splkh} and the \code{splkhfix} functions estimate the odds ratio parameters in the semiparametric odds
#' ratio model using the maximum semiparametric likelihood approach proposed in Chen et al. (2015). The \code{sppenlkh}
#' function performs the network selection using the penalized sepiparametric likelihood approach based on the
#' semiparametric odds ratio model. Currently, only the bilinear log-odds ratio function is allowed.
#'
#' @aliases splkh splkhfix sppenlkh
#'
#' @param y the outcome data matrix. This must be specified.
#' @param x the covariate data matrix. This must be specified.
#' @param fixstruct a matrix indicating which coefficients are set to 0.
#' @param dat the whole data matrix.
#' @param group a vector denoting the sizes of the variables in each group the data matrix divided.
#' @param lambda a vector of penalty parameters. A network is selected for each lambda value.
#' @param niter the maximum number of iterations.
#' @param eps convergence criterion.
#' @param vect the method of vectoring parameter matrix. Either column-wise or row-wise.
#'
#' @usage splkh(y, x, niter = 50, eps = 1e-06, vect = "col")
#' @usage splkhfix(y, x, fixstruct, niter = 50, eps = 1e-06, vect = "col")
#' @usage sppenlkh(dat, group, lambda, niter = 50, eps = 1e-06)
#'
#' @details
#' The \code{splkhfix} function is more general than the \code{splkh} function since the former not only allows for
#' the fixed zero coefficients, it also computes the logarithm of the maximum semiparametric likelihood values.
#'
#' While the outcome data matrix and the variate data matrix have to be specified separately in the \code{splkh} and
#' the \code{splkhfix} functions, for the \code{sppenlkh} function, one has to specify the whole data matrix including
#' the outcome and the covariates.
#'
#' For the \code{splkhfix} function, one has to specify a matrix for a fixed structure, that indicating which
#' coefficients are set to zero. For the \code{sppenlkh} function, a \code{lambda} parameter that specifies the
#' penalties on the semiparametric likelihood, is required.
#'
#' @return
#' The \code{splkh} function returns
#'
#' 1. the odds ratio parameter estimates.
#'
#' 2. the estimated variance-covariance matrix of the odds ratio estimator.
#'
#' The \code{splkhfix} function returns
#'
#' 1. the odds ratio parameter estimates.
#'
#' 2. the estimated variance-covariance matrix of the odds ratio estimator.
#'
#' 3. logarithm of the maximum semiparametric likelihood.
#'
#' The \code{sppenlkh} returns
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
#' # splkh
#' splkh(y = outcomes, x = covariates)
#'
#' # splkhfix
#' structure <- matrix(rbind(p * q), ncol = q)
#' splkhfix(y = outcomes, x = covariates, fixstruct = structure)
#'
#' # sppenlkh
#' n <- 200; p <- 10
#' datmat <- matrix(rnorm(n * p), ncol = p)
#' vargroup <- c(2, 3, 1, 4)
#' penaltyparam <- c(10, 50, 500)
#' sppenlkh(dat = datmat, group = vargroup, lambda = penaltyparam)
#'
#' @references ...
#' @references ...
#'
# Simple approach
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
# Semiparametric likelihood approach with fixed structure
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

# 2. Penalized semiparametric likelihood approach for selecting network connections
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

