#' Multivariate QAR(1) with K = 1
#' 
#' @importFrom stats cov
#' 
#' @description 
#' First series
#' \deqn{\eta_{1,1}(\tau) = F_{a_1,b_1}(\tau)}
#' \deqn{\eta_{2,1}(\tau) = F_{a_2,b_2}(\tau)}
#' Second series
#' \deqn{\eta_{1,2}(\tau) = F_{a_3,b_3}(\tau)}
#' \deqn{\eta_{2,2}(\tau) = F_{a_4,b_4}(\tau)}
#' 
#' @param Y ARRAY, rows are independent, columns have autoregression, 
#'   third dim (length 2) have copula dep (multivariate)
#' @param n.sims,n.thin,n.burnin,n.report (i) Number of iterations after 
#'   burn-in. (ii) Thinning rate. (iii) Number of iterations discarded at the 
#'   beginning. (iv) Report the number of iterations rate.
#' @param inits Initial values
#' @param prior Number. Standard deviation of the normal prior for \eqn{\log a}'s and
#'   \eqn{\log b}'s
#' @param tol Tolerance in the univariate rootfinder. A value too small can lead 
#'   to numerical overflow in the likelihood 
#' @return MATRIX, rows are simulations, columns are parameters
#' \deqn{a_1,b_1,a_2,b_2,a_3,b_3,a_4,b_4,\rho}
#' @export 
MQAR1K1 <- function(Y, 
                    n.sims = 10000,
                    n.thin = 10,
                    n.burnin = 10000,
                    n.report = 1000, 
                    inits = rep(0, 9), 
                    prior = 3, 
                    tol = 10e-12) {
  
  T <- dim(Y)[1]
  L <- dim(Y)[2]
  d <- 9
  sd <- 2.4^2 / d
  epsilon <- 1e-06
  params <- inits
  params <- c(params, logfallMQAR1K1(T, L, inits, Y, prior, tol))
  keepBurnin <- matrix(nrow = 10000 + n.burnin, ncol = d)
  keep   <- matrix(nrow = n.sims / n.thin, ncol = d)
  colnames(keep) <- c("a1", "b1", "a2", "b2", "a3", "b3", "a4", "b4", "rho")
  I <- epsilon * diag(d)
  
  # for 1
  for (b in 1:25) {
    params <- rwBmetropolisMQAR1K1(params, sd * I, T, L, Y, prior, tol)
    
    keepBurnin[b, ] <- params[-10]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[1:25, ]),
    stats::cov(keepBurnin[1:25, ])
  )
  
  # for 2
  for (b in 26:100) {
    params <- rwBmetropolisMQAR1K1(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-10], Sigma[[1]], Sigma[[2]], b - 20)
    
    keepBurnin[b, ] <- params[-10]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[51:100, ]),
    stats::cov(keepBurnin[51:100, ])
  )
  
  # for 3
  for (b in 101:400) {
    params <- rwBmetropolisMQAR1K1(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-10], Sigma[[1]], Sigma[[2]], b - 50)
    
    keepBurnin[b, ] <- params[-10]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[301:400, ]),
    stats::cov(keepBurnin[301:400, ])
  )
  
  # for 4
  for (b in 401:1000) {
    params <- rwBmetropolisMQAR1K1(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-10], Sigma[[1]], Sigma[[2]], b - 300)
    
    keepBurnin[b, ] <- params[-10]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[751:1000, ]),
    stats::cov(keepBurnin[751:1000, ])
  )
  
  # for 5
  for (b in 1001:2000) {
    params <- rwBmetropolisMQAR1K1(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-10], Sigma[[1]], Sigma[[2]], b - 750)
    
    keepBurnin[b, ] <- params[-10]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[1501:2000, ]),
    stats::cov(keepBurnin[1501:2000, ])
  )
  
  # for 6
  for (b in 2001:10000) {
    params <- rwBmetropolisMQAR1K1(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-10], Sigma[[1]], Sigma[[2]], b - 1500)
    
    keepBurnin[b, ] <- params[-10]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[6001:10000, ]),
    stats::cov(keepBurnin[6001:10000, ])
  )
  
  # burn-in
  if (n.burnin < 1000) n.burnin <- 1000
  for (b in 10001:(10000 + n.burnin)) {
    params <- rwBmetropolisMQAR1K1(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-10], Sigma[[1]], Sigma[[2]], b - 9950)
    
    keepBurnin[b, ] <- params[-10]
  }
  
  Sigma <- sd * (stats::cov(keepBurnin[10000 + (1:n.burnin), ]) + I)
  
  print("Burn-in finished")
  
  # FINAL
  for (b in 1:n.sims) {
    if (b %% n.report == 0) print(paste("Iteration:", b))
    
    params <- rwBmetropolisMQAR1K1(params, Sigma, T, L, Y, prior, tol)
    
    # Sigma <- MuSigmaUpdate(params[-10], Sigma[[1]], Sigma[[2]], b + 50)
    
    if (b %% n.thin == 0) {
      keep[b / n.thin, ] <- params[-10]
    }
  }
  
  keep <- exp(keep)
  keep[,9] <- (keep[,9] - 1) / (1 + keep[,9])
  
  return(keep)
}