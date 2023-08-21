#' Koenker and Xiao (2006) QAR(1) model
#' 
#' @importFrom stats cov
#' 
#' @description 
#' This function fits the QAR model by Koenker and Xiao (2006) in a Bayesian 
#' framework, i.e.,
#' \deqn{Q_{Y_{tl}}(\tau \mid y_{t,l-1}) = \mu + \sigma \Phi(\tau) + \min\{\gamma_0 + \gamma_1 \tau, 1\} y_{t,l-1}}
#' with \eqn{\gamma_0 \in (0,1)} and \eqn{\gamma_1 > 0}.
#' 
#' @param Y MATRIX, rows are independent, columns have autoregression
#' @param n.sims,n.thin,n.burnin,n.report (i) Number of iterations after 
#'   burn-in. (ii) Thinning rate. (iii) Number of iterations discarded at the 
#'   beginning. (iv) Report the number of iterations rate.
#' @param inits Initial values (transformed scale \eqn{(-\infty,\infty)})
#' @param prior Vector (length 3). Standard deviation of the zero-mean normal 
#'   prior for (1) \eqn{\mu}, (2) \eqn{\log \sigma}, and (3) 
#'   \eqn{\log \gamma_1}, respectively. 
#'   (\eqn{\gamma_0 \sim U(0,1)} by default.)
#' @param tol Tolerance in the univariate rootfinder. A value too small can 
#'   lead to numerical overflow in the likelihood 
#' @return A \code{"KX2006"} list with elements:
#'   \item{params}{Matrix where rows are simulations and cols are parameters
#'   \deqn{\mu,\sigma,\gamma_0,\gamma_1}}
#'   \item{\code{y}}{Data fitted}
#' @references 
#' Koenker, R. and Xiao, Z. (2006). 
#' “Quantile autoregression.” 
#' \emph{Journal of the American Statistical Association}, \strong{101}(475): 980–990
#' @export 
KX2006 <- function(Y, 
                   n.sims = 10000,
                   n.thin = 10,
                   n.burnin = 10000,
                   n.report = 1000,
                   inits = rep(0, 4), 
                   prior = c(10, 3, 3), 
                   tol = 10e-16) {
  
  T <- nrow(Y)
  L <- ncol(Y)
  d <- 4
  sd <- 2.4^2 / d
  epsilon <- 1e-06
  params <- inits
  params <- c(params, logfallKX2006(T, L, inits, Y, prior, tol))
  keepBurnin <- matrix(nrow = 10000 + n.burnin, ncol = d)
  keep       <- matrix(nrow = n.sims / n.thin, ncol = d)
  colnames(keep) <- c("mu", "sigma", "gamma0", "gamma1")
  I <- epsilon * diag(d)
  
  # for 1
  for (b in 1:25) {
    params <- rwBmetropolisKX2006(params, sd * I, T, L, Y, prior, tol)
    
    keepBurnin[b, ] <- params[-5]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[1:25, ]),
    stats::cov(keepBurnin[1:25, ])
  )
  
  # for 2
  for (b in 26:100) {
    params <- rwBmetropolisKX2006(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-5], Sigma[[1]], Sigma[[2]], b - 20)
    
    keepBurnin[b, ] <- params[-5]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[51:100, ]),
    stats::cov(keepBurnin[51:100, ])
  )
  
  # for 3
  for (b in 101:400) {
    params <- rwBmetropolisKX2006(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-5], Sigma[[1]], Sigma[[2]], b - 50)
    
    keepBurnin[b, ] <- params[-5]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[301:400, ]),
    stats::cov(keepBurnin[301:400, ])
  )
  
  # for 4
  for (b in 401:1000) {
    params <- rwBmetropolisKX2006(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-5], Sigma[[1]], Sigma[[2]], b - 300)
    
    keepBurnin[b, ] <- params[-5]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[751:1000, ]),
    stats::cov(keepBurnin[751:1000, ])
  )
  
  # for 5
  for (b in 1001:2000) {
    params <- rwBmetropolisKX2006(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-5], Sigma[[1]], Sigma[[2]], b - 750)
    
    keepBurnin[b, ] <- params[-5]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[1501:2000, ]),
    stats::cov(keepBurnin[1501:2000, ])
  )
  
  # for 6
  for (b in 2001:10000) {
    params <- rwBmetropolisKX2006(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
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
    params <- rwBmetropolisKX2006(params, sd * (Sigma[[2]] + I), T, L, Y, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-5], Sigma[[1]], Sigma[[2]], b - 9950)
    
    keepBurnin[b, ] <- params[-5]
  }
  
  Sigma <- sd * (stats::cov(keepBurnin[10000 + (1:n.burnin), ]) + I)
  
  print("Burn-in finished")
  
  # FINAL
  for (b in 1:n.sims) {
    if (b %% n.report == 0) print(paste("Iteration:", b))
    
    params <- rwBmetropolisKX2006(params, Sigma, T, L, Y, prior, tol)
    
    # Sigma <- MuSigmaUpdate(params[-5], Sigma[[1]], Sigma[[2]], b + 50)
    
    if (b %% n.thin == 0) {
      keep[b / n.thin, ] <- params[-5]
    }
  }
  
  keep[,2:4] <- exp(keep[,2:4]) 
  keep[,3] <- keep[,3] / (1 + keep[,3])
  
  keep.list <- list()
  keep.list$params <- keep
  keep.list$y <- Y

  class(keep.list) <- "KX2006"

  return(keep.list)
}