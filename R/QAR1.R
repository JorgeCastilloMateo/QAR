#' QAR(1) with K = 1
#' 
#' @importFrom stats cov
#' 
#' @description 
#' This function fits the model 
#' \deqn{Q_{Y_{tl}}(\tau \mid y_{t,l-1}) = y_{t,l-1} \eta_1(\tau) + (1 - y_{t,l-1}) \eta_2(\tau)}
#' with
#' \deqn{\eta_1(\tau) = F(\tau \mid a_1,b_1)}
#' \deqn{\eta_2(\tau) = F(\tau \mid a_2,b_2)}
#' and \eqn{F(\tau \mid a,b)} the Kumaraswamy cdf with parameters \eqn{a} and \eqn{b}
#' 
#' @param Y MATRIX, rows are independent, columns have autoregression
#' @param n.sims,n.thin,n.burnin,n.report (i) Number of iterations after 
#'   burn-in. (ii) Thinning rate. (iii) Number of iterations discarded at the 
#'   beginning. (iv) Report the number of iterations rate.
#' @param inits Initial values (transformed scale \eqn{(-\infty,\infty)})
#' @param prior Number. Standard deviation of the zero-mean normal prior for 
#'   \eqn{\log a}'s and \eqn{\log b}'s
#' @param tol Tolerance in the univariate rootfinder. A value too small can lead 
#'   to numerical overflow in the likelihood 
#' @return A \code{"QAR1K1"} list with elements:
#'   \item{params}{Matrix where rows are simulations and cols are parameters
#'   \deqn{a_1,b_1,a_2,b_2}}
#'   \item{\code{y}}{Data fitted}
#' @export 
QAR1K1 <- function(Y, 
                   n.sims = 10000,
                   n.thin = 10,
                   n.burnin = 10000,
                   n.report = 1000, 
                   inits = rep(0, 4),
                   prior = 3, 
                   tol = 10e-16) {
  
  T <- nrow(Y)
  L <- ncol(Y)
  d <- 4
  sd <- 2.4^2 / d
  epsilon <- 1e-06
  params <- inits
  params <- c(params, logfallQAR1K1(T, L, inits, Y, prior, tol))
  keepBurnin <- matrix(nrow = 10000 + n.burnin, ncol = d)
  keep       <- matrix(nrow = n.sims / n.thin, ncol = d)
  colnames(keep) <- c("a1", "b1", "a2", "b2")
  I <- epsilon * diag(d)
  
  # for 1
  for (b in 1:25) {
    params <- rwBmetropolisQAR1K1(params, sd * I, T, L, Y, prior, tol)
    
    keepBurnin[b, ] <- params[-5]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[1:25, ]),
    stats::cov(keepBurnin[1:25, ])
  )
  
  # for 2
  for (b in 26:100) {
    params <- rwBmetropolisQAR1K1(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-5], Sigma[[1]], Sigma[[2]], b - 20)
    
    keepBurnin[b, ] <- params[-5]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[51:100, ]),
    stats::cov(keepBurnin[51:100, ])
  )
  
  # for 3
  for (b in 101:400) {
    params <- rwBmetropolisQAR1K1(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-5], Sigma[[1]], Sigma[[2]], b - 50)
    
    keepBurnin[b, ] <- params[-5]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[301:400, ]),
    stats::cov(keepBurnin[301:400, ])
  )
  
  # for 4
  for (b in 401:1000) {
    params <- rwBmetropolisQAR1K1(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-5], Sigma[[1]], Sigma[[2]], b - 300)
    
    keepBurnin[b, ] <- params[-5]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[751:1000, ]),
    stats::cov(keepBurnin[751:1000, ])
  )
  
  # for 5
  for (b in 1001:2000) {
    params <- rwBmetropolisQAR1K1(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-5], Sigma[[1]], Sigma[[2]], b - 750)
    
    keepBurnin[b, ] <- params[-5]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[1501:2000, ]),
    stats::cov(keepBurnin[1501:2000, ])
  )
  
  # for 6
  for (b in 2001:10000) {
    params <- rwBmetropolisQAR1K1(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-5], Sigma[[1]], Sigma[[2]], b - 1500)
    
    keepBurnin[b, ] <- params[-5]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[6001:10000, ]),
    stats::cov(keepBurnin[6001:10000, ])
  )
  
  # burn-in
  if (n.burnin < 1000) n.burnin <- 1000
  for (b in 10001:(10000 + n.burnin)) {
    params <- rwBmetropolisQAR1K1(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-5], Sigma[[1]], Sigma[[2]], b - 9950)
    
    keepBurnin[b, ] <- params[-5]
  }
  
  Sigma <- sd * (stats::cov(keepBurnin[10000 + (1:n.burnin), ]) + I)
  
  print("Burn-in finished")
  
  # FINAL
  for (b in 1:n.sims) {
    if (b %% n.report == 0) print(paste("Iteration:", b))
    
    params <- rwBmetropolisQAR1K1(params, Sigma, T, L, Y, prior, tol)
    
    #Sigma <- MuSigmaUpdate(params[-5], Sigma[[1]], Sigma[[2]], b + 50)
    
    if (b %% n.thin == 0) {
      keep[b / n.thin, ] <- params[-5]
    }
  }
  
  keep <- exp(keep)
  
  keep.list <- list()
  keep.list$params <- keep
  keep.list$y <- Y
  
  class(keep.list) <- "QAR1K1"
  
  return(keep.list)
}

#' QAR(1) with K = 2
#' 
#' @importFrom stats cov
#' 
#' @description 
#' This function fits the model 
#' \deqn{Q_{Y_{tl}}(\tau \mid y_{t,l-1}) = y_{t,l-1} \eta_1(\tau) + (1 - y_{t,l-1}) \eta_2(\tau)}
#' with
#' \deqn{\eta_1(\tau) = \lambda_1 F(\tau \mid a_1,b_1) + (1 - \lambda_1) F(\tau \mid a_2,b_2)}
#' \deqn{\eta_2(\tau) = \lambda_2 F(\tau \mid a_3,b_3) + (1 - \lambda_2) F(\tau \mid a_4,b_4)}
#' and \eqn{F(\tau \mid a,b)} the Kumaraswamy cdf with parameters \eqn{a} and \eqn{b}
#' 
#' @param Y MATRIX, rows are independent, columns have autoregression
#' @param n.sims,n.thin,n.burnin,n.report (i) Number of iterations after 
#'   burn-in. (ii) Thinning rate. (iii) Number of iterations discarded at the 
#'   beginning. (iv) Report the number of iterations rate.
#' @param inits Initial values (transformed scale \eqn{(-\infty,\infty)})
#' @param prior Number. Standard deviation of the zero-mean normal prior for 
#'   \eqn{\log a}'s and \eqn{\log b}'s. 
#'   (\eqn{\lambda_1,\lambda_2 \sim U(0,1/2)} by default.)
#' @param tol Tolerance in the univariate rootfinder. A value too small can 
#'   lead to numerical overflow in the likelihood 
#' @return A \code{"QAR1K2"} list with elements:
#'   \item{params}{Matrix where rows are simulations and cols are parameters
#'   \deqn{a_1,b_1,a_2,b_2,a_3,b_3,a_4,b_4,\lambda_1,\lambda_2}}
#'   \item{\code{y}}{Data fitted}
#' @export 
QAR1K2 <- function(Y, 
                   n.sims = 10000,
                   n.thin = 10,
                   n.burnin = 10000,
                   n.report = 1000,
                   inits = rep(0, 10), 
                   prior = 3 / 2, 
                   tol = 10e-16) {
  
  T <- nrow(Y)
  L <- ncol(Y)
  d <- 10
  sd <- 2.4^2 / d
  epsilon <- 1e-06
  params <- inits
  params <- c(params, logfallQAR1K2(T, L, inits, Y, prior, tol))
  keepBurnin <- matrix(nrow = 10000 + n.burnin, ncol = d)
  keep       <- matrix(nrow = n.sims / n.thin, ncol = d)
  colnames(keep) <- c("a1", "b1", "a2", "b2", "a3", "b3", "a4", "b4", "lambda1", "lambda2")
  I <- epsilon * diag(d)
  
  # for 1
  for (b in 1:25) {
    params <- rwBmetropolisQAR1K2(params, sd * I, T, L, Y, prior, tol)
    
    keepBurnin[b, ] <- params[-11]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[1:25, ]),
    stats::cov(keepBurnin[1:25, ])
  )
  
  # for 2
  for (b in 26:100) {
    params <- rwBmetropolisQAR1K2(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-11], Sigma[[1]], Sigma[[2]], b - 20)
    
    keepBurnin[b, ] <- params[-11]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[51:100, ]),
    stats::cov(keepBurnin[51:100, ])
  )
  
  # for 3
  for (b in 101:400) {
    params <- rwBmetropolisQAR1K2(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-11], Sigma[[1]], Sigma[[2]], b - 50)
    
    keepBurnin[b, ] <- params[-11]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[301:400, ]),
    stats::cov(keepBurnin[301:400, ])
  )
  
  # for 4
  for (b in 401:1000) {
    params <- rwBmetropolisQAR1K2(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-11], Sigma[[1]], Sigma[[2]], b - 300)
    
    keepBurnin[b, ] <- params[-11]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[751:1000, ]),
    stats::cov(keepBurnin[751:1000, ])
  )
  
  # for 5
  for (b in 1001:2000) {
    params <- rwBmetropolisQAR1K2(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-11], Sigma[[1]], Sigma[[2]], b - 750)
    
    keepBurnin[b, ] <- params[-11]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[1501:2000, ]),
    stats::cov(keepBurnin[1501:2000, ])
  )
  
  # for 6
  for (b in 2001:10000) {
    params <- rwBmetropolisQAR1K2(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-11], Sigma[[1]], Sigma[[2]], b - 1500)
    
    keepBurnin[b, ] <- params[-11]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[6001:10000, ]),
    stats::cov(keepBurnin[6001:10000, ])
  )
  
  # burn-in
  if (n.burnin < 1000) n.burnin <- 1000
  for (b in 10001:(10000 + n.burnin)) {
    params <- rwBmetropolisQAR1K2(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-11], Sigma[[1]], Sigma[[2]], b - 9950)
    
    keepBurnin[b, ] <- params[-11]
  }
  
  Sigma <- sd * (stats::cov(keepBurnin[10000 + (1:n.burnin), ]) + I)
  
  print("Burn-in finished")
  
  # FINAL
  for (b in 1:n.sims) {
    if (b %% n.report == 0) print(paste("Iteration:", b))
    
    params <- rwBmetropolisQAR1K2(params, Sigma, T, L, Y, prior, tol)
    
    #Sigma <- MuSigmaUpdate(params[-11], Sigma[[1]], Sigma[[2]], b + 50)
    
    if (b %% n.thin == 0) {
      keep[b / n.thin, ] <- params[-11]
    }
  }

  keep <- data.frame(exp(keep))
  keep[,9:10] <- (0.5 * keep[,9:10]) / (1 + keep[,9:10])
  
  keep.list <- list()
  keep.list$params <- keep
  keep.list$y <- Y
  
  class(keep.list) <- "QAR1K2"
  
  return(keep.list)
  
  return(keep)
}

#' QAR(1) with K in a Basis (fixed a's and b's)
#' 
#' @importFrom stats cov
#' 
#' @description 
#' This function fits the model 
#' \deqn{Q_{Y_{tl}}(\tau \mid y_{t,l-1}) = y_{t,l-1} \eta_1(\tau) + (1 - y_{t,l-1}) \eta_2(\tau)}
#' where for \eqn{j=1,2},
#' \deqn{\eta_j(\tau) = \sum_{k=1}^{K} \lambda_k^{(j)} F(\tau \mid a_k,b_k)}
#' with \eqn{\lambda_K^{(j)} = 1 - \sum_{k=1}^{K-1} \lambda_k^{(j)}},
#' \eqn{(a_k,b_k)_{k=1}^{K}} fixed in the arguments, and \eqn{F(\tau \mid a,b)}
#' the Kumaraswamy cdf with parameters \eqn{a} and \eqn{b}
#' 
#' @param Y MATRIX, rows are independent, columns have autoregression
#' @param K integer with the number of terms in the basis
#' @param coef VECTOR of length \eqn{2 \times K} with 
#'   \eqn{(a_1,b_1,\ldots,a_K,b_K)}
#' @param n.sims,n.thin,n.burnin,n.report (i) Number of iterations after 
#'   burn-in. (ii) Thinning rate. (iii) Number of iterations discarded at the 
#'   beginning. (iv) Report the number of iterations rate.
#' @param inits Initial values (transformed scale \eqn{(-\infty,\infty)})
#' @param prior Number. Standard deviation of the zero-mean normal prior for 
#'   the \eqn{\phi}'s 
#'   (the logistic transformation of the \eqn{\lambda}'s)
#' @param tol Tolerance in the univariate rootfinder. A value too small can lead 
#'   to numerical overflow in the likelihood 
#' @return A \code{"QAR1K1"} list with elements:
#'   \item{params}{Matrix where rows are simulations and cols are parameters
#'   \deqn{\lambda_1^{(1)},\ldots,\lambda_{K-1}^{(1)},
#'         \lambda_1^{(2)},\ldots,\lambda_{K-1}^{(2)}}}
#'   \item{\code{y}}{Data fitted}
#'   \item{\code{coef}}{Coefficients \eqn{(a_1,b_1,\ldots,a_K,b_K)}}
#' @export 
QAR1K <- function(Y, 
                  K = 8, 
                  coef = c(.5,.5,1,1,.5,2,2,.5,1,16,16,1,4,8,8,4), 
                  n.sims = 10000,
                  n.thin = 10,
                  n.burnin = 10000,
                  n.report = 1000,
                  inits = rep(0, 14), 
                  prior = 3, 
                  tol = 10e-16) {
  
  if (length(coef) != 2 * K) stop("'length(coef)' must be '2 * K'")
  if (length(inits) != 2 * K - 2) stop("'length(inits)' must be '2 * K - 2'")

  T <- nrow(Y)
  L <- ncol(Y)
  d <- 2 * K - 2
  d1 <- d + 1
  sd <- 2.4^2 / d
  epsilon <- 1e-06
  params <- inits
  params <- c(params, logfallQAR1K(T, L, inits, Y, K, coef, prior, tol))
  keepBurnin <- matrix(nrow = 10000 + n.burnin, ncol = d)
  keep   <- matrix(nrow = n.sims / n.thin, ncol = d)
  colnames(keep) <- c(paste0("l", 1:(K-1), "1"), paste0("l", 1:(K-1), "2"))
  I <- epsilon * diag(d)
  
  # for 1
  for (b in 1:25) {
    params <- rwBmetropolisQAR1K(params, sd * I, T, L, Y, K, coef, prior, tol)
    
    keepBurnin[b, ] <- params[-d1]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[1:25, ]),
    stats::cov(keepBurnin[1:25, ])
  )

  # for 2
  for (b in 26:100) {
    params <- rwBmetropolisQAR1K(params, sd * (Sigma[[2]] + I), T, L, Y, K, coef, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-d1], Sigma[[1]], Sigma[[2]], b - 20)
    
    keepBurnin[b, ] <- params[-d1]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[51:100, ]),
    stats::cov(keepBurnin[51:100, ])
  )
  
  # for 3
  for (b in 101:400) {
    params <- rwBmetropolisQAR1K(params, sd * (Sigma[[2]] + I), T, L, Y, K, coef, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-d1], Sigma[[1]], Sigma[[2]], b - 50)
    
    keepBurnin[b, ] <- params[-d1]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[301:400, ]),
    stats::cov(keepBurnin[301:400, ])
  )
  
  # for 4
  for (b in 401:1000) {
    params <- rwBmetropolisQAR1K(params, sd * (Sigma[[2]] + I), T, L, Y, K, coef, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-d1], Sigma[[1]], Sigma[[2]], b - 300)
    
    keepBurnin[b, ] <- params[-d1]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[751:1000, ]),
    stats::cov(keepBurnin[751:1000, ])
  )
  
  # for 5
  for (b in 1001:2000) {
    params <- rwBmetropolisQAR1K(params, sd * (Sigma[[2]] + I), T, L, Y, K, coef, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-d1], Sigma[[1]], Sigma[[2]], b - 750)
    
    keepBurnin[b, ] <- params[-d1]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[1501:2000, ]),
    stats::cov(keepBurnin[1501:2000, ])
  )
  
  # for 6
  for (b in 2001:10000) {
    params <- rwBmetropolisQAR1K(params, sd * (Sigma[[2]] + I), T, L, Y, K, coef, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-d1], Sigma[[1]], Sigma[[2]], b - 1500)
    
    keepBurnin[b, ] <- params[-d1]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[6001:10000, ]),
    stats::cov(keepBurnin[6001:10000, ])
  )
  
  # burn-in
  if (n.burnin < 1000) n.burnin <- 1000
  for (b in 10001:(10000 + n.burnin)) {
    params <- rwBmetropolisQAR1K(params, sd * (Sigma[[2]] + I), T, L, Y, K, coef, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-d1], Sigma[[1]], Sigma[[2]], b - 9950)
    
    keepBurnin[b, ] <- params[-d1]
  }
  
  Sigma <- sd * (stats::cov(keepBurnin[10000 + (1:n.burnin), ]) + I)
  
  print("Burn-in finished")
  
  # FINAL
  for (b in 1:n.sims) {
    if (b %% n.report == 0) print(paste("Iteration:", b))
    
    params <- rwBmetropolisQAR1K(params, Sigma, T, L, Y, K, coef, prior, tol)
    
    # Sigma <- MuSigmaUpdate(params[-d1], Sigma[[1]], Sigma[[2]], b + 50)
    
    if (b %% n.thin == 0) {
      keep[b / n.thin, ] <- params[-d1]
    }
  }
  
  keep <- exp(keep)
  
  sum1 <- 1 + rowSums(keep[,1:(K-1)])
  sum2 <- 1 + rowSums(keep[,K:(2*K-2)])
  
  keep[,1:(K-1)]   <- keep[,1:(K-1)]   / sum1
  keep[,K:(2*K-2)] <- keep[,K:(2*K-2)] / sum2
  
  keep.list <- list()
  keep.list$params <- keep
  keep.list$y <- Y
  keep.list$coef <- coef
  
  class(keep.list) <- "QAR1K"
  
  return(keep.list)
}
