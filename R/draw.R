#' @importFrom graphics lines
NULL
# for internal use
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
