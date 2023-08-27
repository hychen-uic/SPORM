#' @importFrom graphics lines
NULL

#' This function is used for diagnostic of convergence
#'
#' This function is used to draw the progression of the parameter estimates
#' in the iterations of the permutation likelihood approach to monitor the convergence
#' of the parameter estimates.
#'
#' @param x input argument.
#'
#' @export
#'
#----------------------
draw <- function(x){
  n <- dim(x)[1]
  p <- dim(x)[2]
  low <- min(x) * (1 + 0.2)
  upp=max(x) * (1 + 0.2)
  plot(x = c(-35), y = c(45), xlim = c(0,n+5), ylim = c(low, upp),
       xlab = "Number of Iterations",ylab="Odds Ratio Estimates")
  for(i in 2:p) {for(k in 1:(i-1)) {
    lines(c(1:n), x[1:n, i, k], type = "l", col = 'red')
  }}
}
