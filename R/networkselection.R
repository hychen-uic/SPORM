#' The network selection from a set of selected networks using output from one of
#' the penalized likelihood approaches
#'
#' This approach selects the network using the BIC.
#'
#' @param dat a data matrix of n*node dimension
#' @param group a vector of positive integers of length ng such that sum(group) = node
#' @param lambda a vector of penalty values for network selection, each penalty value determines a network
#' @param network a set of networks determined by one of the penalized approaches
#' @param method the method used in the likelihood approach ("pw", "sp", or "pm")
#' @param criterion selection criterion such as BIC
#'
#' @details This approach selects the network using the BIC, by refitting each
#' network to obtain the log-likelihood (objective function) and then
#' uses BIC for the selection.
#'
#' @return Estimate of the model parameters (vectorized) and
#'         the covariance matrix estimate corresponding to the vectorized parameters.
#'
#' @references Chen, H.Y. (2022). Semiparametric Odds Ratio Model and its Application. CRC press.
#' @references Chen, H. Y. and Chen, J. (2020). Network selection through semiparametric odds ratio model. Manuscript.
#'
#' @examples \dontrun{
#' n <- 200; p <- 10
#' datmat <- matrix(rnorm(n * p), ncol = p)
#' vargroup <- c(2, 3, 1, 4)
#' penaltyparam <- c(10, 50, 500)
#' fit <- pmpenlkh(dat = datmat, group = vargroup, lambda = penaltyparam,
#'                 nsamp = 2e4, niter = 50)
#' network <- fit[[4]]
#' fit0 <- networkselect(dat = datmat, group = vargroup, lambda = penaltyparam,
#'                       network = network, method = 'pm', criterion = 'BIC')
#' fit0[[2]]
#' }
#'
#' @export
#'
networkselect <- function(dat, group, lambda, network, method, criterion) {
  n <- dim(dat)[1]
  p <- dim(dat)[2]
  ng <- length(group)
  nlam <- length(lambda)
  q <- p * (p - 1)/2
  for (i in 1:ng){
    q <- q - group[i] * (group[i] - 1)/2
  }
  if(nlam != dim(network)[3]) {
    # change exit to stop
    stop("inconsistent length of lambda and the third dimension of network")
    # exit
  } else{
    loglkh <- rep(0,nlam)
    for(i in 1:nlam) {
      if(method == 'pw') {
        kk <- 0
        for(k in 1:ng) {

          if(k == 1) {
            y <- dat[, 1:group[k]]
            x <- dat[, (group[k] + 1):p]
            struct <- network[1:group[k], (group[k]+1):p, i]
            loglkh[i] <- pwlkhfix(y = y,x = x,fixstruct = struct)[[3]]
          } else{if(k < ng) {
            y <- dat[, (kk + 1):(kk + group[k])]
            x <- cbind(dat[,1:kk], dat[, (kk + group[k] + 1):p])
            struct <- cbind(t(network[(kk+1):(kk+group[k]), 1:kk, i]),
                            network[1:kk, (kk + group[k] + 1):p, i])
            loglkh[i]<- loglkh[i] + pwlkhfix(y = y, x = x, fixstruct = struct)[[3]]
          } else {
            y <- dat[, (p - group[ng] + 1):p]
            x <- dat[, 1:(p - group[ng])]
            struct <- t(network[1:(p - group[ng]), (p - group[ng] + 1):p, i])
            loglkh[i] <- loglkh[i] + pwlkhfix(y = y, x = x, fixstruct = struct)[[3]]
          }
          }
          kk <- kk + group[k]
        }
        loglkh[i] <- loglkh[i]/2 #for repeatedly counting
      }
      if(method == 'sp'){
        kk <- 0
        for(k in 1:ng){
          if(k == 1){
            y <- dat[, 1:group[k]]
            x <- dat[, (group[k] + 1):p]
            struct <- network[1:group[k], (group[k]+1):p, i]
            loglkh[i] <- splkhfix(y = y, x = x, fixstruct = struct)[[3]]
          } else{if(k < ng) {
            y <- dat[, (kk + 1):(kk + group[k])]
            x <- cbind(dat[, 1:kk], dat[, (kk + group[k] + 1):p])
            struct <- cbind(t(network[(kk + 1):(kk + group[k]), 1:kk, i]),
                            network[1:kk, (kk + group[k] + 1):p, i])
            loglkh[i] <- loglkh[i] + splkhfix(y = y, x = x, fixstruct = struct)[[3]]
          } else{
            y <- dat[, (p - group[ng] + 1):p]
            x <- dat[, 1:(p - group[ng])]
            struct <- t(network[1:(p - group[ng]), (p - group[ng] + 1):p, i])
            loglkh[i] <- loglkh[i] + splkhfix(y = y, x = x, fixstruct = struct)[[3]]
          }
          }
          kk <- kk + group[k]
        }
        loglkh[i] <- loglkh[i]/2 #for repeatedly counting
      }
      if(method == 'pm'){
        loglkh[i] <- pmlkhfix(dat=dat,group=group,fixstruct=network[,,i])[[3]]
      }
    }}

  if(criterion == 'BIC'){
    BIC <- rep(0, nlam)

    selmodel <- 1
    BIC[1] <- -2 * loglkh[1] + (sum(network[ , , 1])/2) * log(n)
    minBIC <- BIC[1]
    for(i in 2:nlam){
      BIC[i] <- -2 * loglkh[i] + (sum(network[, , i])/2) * log(n)
      if(BIC[i]<minBIC){
        minBIC <- BIC[i]
        selmodel <- i
      }
    }
  }
  return(list(selmodel, network[ , , selmodel], BIC, loglkh))
}
