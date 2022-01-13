#' Network selection
#'
#' This function selects one network among a set of networks using a specified criterion based on a particular
#' fitting approach.
#'
#' @param dat the whole data matrix.
#' @param group a vector indicating the sizes of the variables in each group the data matrix divides.
#' @param lambda a vector of penalty parameters. A network is selected for each lambda value.
#' @param network a set of networks correspond to lambda vector.
#' @param method the method used to fit the network model, choices include: "pw" for pairwise pseudo-likelihood
#' approach, "pm" for permutation likelihood approach, and "sp" for semiparametric likelihood approach.
#' @param criterion the network selection criterion: "BIC".
#'
#' @details
#' One has to specify the whole data matrix including the outcome and covariates, a vector indicating the sizes of
#' variables in each group that data matrix divides has to be specified.
#'
#' One can choose the desired approach, which including "pw", "pm", or "sp" for fitting the network model. Alsom a network
#' selection criterion has to be specified. Currently, the implemented network selection criterion is "BIC".
#'
#' @return
#' The function returns
#'
#' 1. the index of the selected model
#'
#' 2. the structure of the selected network
#'
#' 3. all the BICs corresponding to the set of networks
#'
#' 4. all the log-likelihood values corresponding to the set of networks: sum of the node-wise pseudo-likelihoods for
#' \code{method = "pw"}; sum of the node-wise semiparametric likelihoods for \code{method = "sp"}; the joint permutation
#' likelihoods for \code{method = "pm"}
#'
#' @examples
#' \dontrun{
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
#' @references Chen, H. Y. and Chen, J. (2021). Semiparametric odds ratio model for network detection. Manuscript.
#' @references Chen, H. Y. (2021). Semiparametric odds ratio model and its applications. draft.
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
