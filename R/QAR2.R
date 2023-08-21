#' QAR(2) with K = 1
#' 
#' @importFrom stats cov
#' 
#' @description 
#' This function fits the model 
#' \deqn{Q_{Y_{tl}}(\tau \mid y_{t,l-1}, y_{t,l-2}) = \pi y_{t,l-1} \eta_1(\tau) + (1 - \pi) y_{t,l-2} \eta_2(\tau) + (1 - \pi y_{t,l-1} - (1 - \pi) y_{t,l-2}) \eta_3(\tau)}
#' with
#' \deqn{\eta_1(\tau) = F(\tau \mid a_1,b_1)}
#' \deqn{\eta_2(\tau) = F(\tau \mid a_2,b_2)}
#' \deqn{\eta_3(\tau) = F(\tau \mid a_3,b_3)}
#' and \eqn{F(\tau \mid a,b)} the Kumaraswamy cdf with parameters \eqn{a} and \eqn{b}
#' 
#' @param Y MATRIX, rows are independent, columns have autoregression
#' @param n.sims,n.thin,n.burnin,n.report (i) Number of iterations after 
#'   burn-in. (ii) Thinning rate. (iii) Number of iterations discarded at the 
#'   beginning. (iv) Report the number of iterations rate.
#' @param inits Initial values (transformed scale \eqn{(-\infty,\infty)})
#' @param prior Number. Standard deviation of the zero-mean  normal prior for 
#'   \eqn{\log a}'s and \eqn{\log b}'s. 
#'   (\eqn{\pi \sim U(0,1)} by default.)
#' @param tol Tolerance in the univariate rootfinder. A value too small can 
#'   lead to numerical overflow in the likelihood 
#' @return A \code{"QAR2K1"} list with elements:
#'   \item{params}{Matrix where rows are simulations and cols are parameters
#'   \deqn{a_1,b_1,a_2,b_2,a_3,b_3,\pi}}
#'   \item{\code{y}}{Data fitted}
#' @export 
QAR2K1 <- function(Y, 
                   n.sims = 10000,
                   n.thin = 10,
                   n.burnin = 10000,
                   n.report = 1000, 
                   inits = rep(0, 7), 
                   prior = 3 / 2, 
                   tol = 10e-16) {
  
  T <- nrow(Y)
  L <- ncol(Y)
  d <- 7
  sd <- 2.4^2 / d
  epsilon <- 1e-06
  params <- inits
  params <- c(params, logfallQAR2K1(T, L, inits, Y, prior, tol))
  keepBurnin <- matrix(nrow = 10000 + n.burnin, ncol = d)
  keep       <- matrix(nrow = n.sims / n.thin, ncol = d)
  colnames(keep) <- c("a1", "b1", "a2", "b2", "a3", "b3", "pi")
  I <- epsilon * diag(d)
  
  # for 1
  for (b in 1:25) {
    params <- rwBmetropolisQAR2K1(params, sd * I, T, L, Y, prior, tol)
    
    keepBurnin[b, ] <- params[-8]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[1:25, ]),
    stats::cov(keepBurnin[1:25, ])
  )
  
  # for 2
  for (b in 26:100) {
    params <- rwBmetropolisQAR2K1(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-8], Sigma[[1]], Sigma[[2]], b - 20)
    
    keepBurnin[b, ] <- params[-8]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[51:100, ]),
    stats::cov(keepBurnin[51:100, ])
  )
  
  # for 3
  for (b in 101:400) {
    params <- rwBmetropolisQAR2K1(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-8], Sigma[[1]], Sigma[[2]], b - 50)
    
    keepBurnin[b, ] <- params[-8]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[301:400, ]),
    stats::cov(keepBurnin[301:400, ])
  )
  
  # for 4
  for (b in 401:1000) {
    params <- rwBmetropolisQAR2K1(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-8], Sigma[[1]], Sigma[[2]], b - 300)
    
    keepBurnin[b, ] <- params[-8]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[751:1000, ]),
    stats::cov(keepBurnin[751:1000, ])
  )
  
  # for 5
  for (b in 1001:2000) {
    params <- rwBmetropolisQAR2K1(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-8], Sigma[[1]], Sigma[[2]], b - 750)
    
    keepBurnin[b, ] <- params[-8]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[1501:2000, ]),
    stats::cov(keepBurnin[1501:2000, ])
  )
  
  # for 6
  for (b in 2001:10000) {
    params <- rwBmetropolisQAR2K1(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-8], Sigma[[1]], Sigma[[2]], b - 1500)
    
    keepBurnin[b, ] <- params[-8]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[6001:10000, ]),
    stats::cov(keepBurnin[6001:10000, ])
  )
  
  # burn-in
  if (n.burnin < 1000) n.burnin <- 1000
  for (b in 10001:(10000 + n.burnin)) {
    params <- rwBmetropolisQAR2K1(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-8], Sigma[[1]], Sigma[[2]], b - 9950)
    
    keepBurnin[b, ] <- params[-8]
  }
  
  Sigma <- sd * (stats::cov(keepBurnin[10000 + (1:n.burnin), ]) + I)
  
  print("Burn-in finished")
  
  # FINAL
  for (b in 1:n.sims) {
    if (b %% n.report == 0) print(paste("Iteration:", b))
    
    params <- rwBmetropolisQAR2K1(params, Sigma, T, L, Y, prior, tol)
    
    # Sigma <- MuSigmaUpdate(params[-8], Sigma[[1]], Sigma[[2]], b + 50)
    
    if (b %% n.thin == 0) {
      keep[b / n.thin, ] <- params[-8]
    }
  }
  
  keep <- exp(keep)
  keep[,7] <- keep[,7] / (1 + keep[,7])
  
  keep.list <- list()
  keep.list$params <- keep
  keep.list$y <- Y
  
  class(keep.list) <- "QAR2K1"
  
  return(keep.list)
}