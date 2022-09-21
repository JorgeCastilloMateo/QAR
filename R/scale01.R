#' Scale in \eqn{(0,1)}
#' 
#' @description This function scales a variable in \eqn{(0,1)} following ideas from the 
#' theory of order statistics in uniform random variables
#' \deqn{y_{t} = \frac{y_t^* - m}{M - m},}
#' where \eqn{m = (T y_{(1)}^{*} - y_{(T)}^{*})/(T-1)} and 
#' \eqn{M = (T y_{(T)}^{*} - y_{(1)}^{*})/(T-1)}.
#' 
#' @param y VECTOR or MATRIX to scale between 0 and 1
#' @return VECTOR or MATRIX with values between 0 and 1
#' @export 
scale01 <- function(y) {
  T <- length(y)
  minY <- min(y)
  maxY <- max(y)
  m <- (T * minY - maxY) / (T - 1)
  M <- (T * maxY - minY) / (T - 1)
  y <- (y - m) / (M - m)
  return(y)
}