#' Spatial QAR(1) with K = 1
#' 
#' @importFrom stats cov 
#' @importFrom stats dist
#' 
#' @description 
#' \deqn{\eta_1(\tau) = F_{a_1,b_1}(\tau)}
#' \deqn{\eta_2(\tau) = F_{a_2,b_2}(\tau)}
#' 
#' Fast convergence
#' 
#' @param Y ARRAY, rows are independent, columns have autoregression, 
#'   third dim have copula dep (space)
#' @param coords COORDINATES \eqn{n \times 2} of the locations
#' @param n.sims,n.thin,n.burnin,n.report (i) Number of iterations after 
#'   burn-in. (ii) Thinning rate. (iii) Number of iterations discarded at the 
#'   beginning. (iv) Report the number of iterations rate.
#' @param inits Initial values
#' @param prior Number. Standard deviation of the normal prior for \eqn{\log a}'s and
#'   \eqn{\log b}'s
#' @param tol Tolerance in the univariate rootfinder. A value too small can lead 
#'   to numerical overflow in the likelihood 
#' @return MATRIX, rows are simulations, columns are parameters
#' \deqn{a_1,b_1,a_2,b_2,\gamma}
#' @export 
SQAR1K1 <- function(Y, 
                    coords, 
                    n.sims = 10000,
                    n.thin = 10,
                    n.burnin = 10000,
                    n.report = 1000, 
                    inits = rep(0, 5), 
                    prior = 3, 
                    tol = 10e-16) {
  
  T <- dim(Y)[1]
  L <- dim(Y)[2]
  n <- dim(Y)[3]
  d <- 5
  sd <- 2.4^2 / d
  epsilon <- 1e-06
  params <- inits
  keepBurnin <- matrix(nrow = 10000 + n.burnin, ncol = d)
  keep       <- matrix(nrow = n.sims / n.thin, ncol = d)
  colnames(keep) <- c("a1", "b1", "a2", "b2", "gamma")
  I <- epsilon * diag(d)
  
  distance <- stats::dist(coords)
  phi <- 3 / max(distance) # fixed phi
  dAux <- matrix(0, nrow = n, ncol = n)
  dAux[lower.tri(dAux)] <- distance
  dAux <- dAux + t(dAux)
  R.phi <- exp(- phi * dAux)
  Id <- diag(n)
  
  params <- c(params, logfallSQAR1K1(T, L, n, inits, Y, R.phi, Id, prior, tol))
  
  # for 1
  for (b in 1:25) {
    params <- rwBmetropolisSQAR1K1(params, sd * I, T, L, n, Y, R.phi, Id, prior, tol)
    
    keepBurnin[b, ] <- params[-6]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[1:25, ]),
    stats::cov(keepBurnin[1:25, ])
  )
  
  # for 2
  for (b in 26:100) {
    params <- rwBmetropolisSQAR1K1(params, sd * (Sigma[[2]] + I), T, L, n, Y, R.phi, Id, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-6], Sigma[[1]], Sigma[[2]], b - 20)
    
    keepBurnin[b, ] <- params[-6]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[51:100, ]),
    stats::cov(keepBurnin[51:100, ])
  )
  
  # for 3
  for (b in 101:400) {
    params <- rwBmetropolisSQAR1K1(params, sd * (Sigma[[2]] + I), T, L, n, Y, R.phi, Id, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-6], Sigma[[1]], Sigma[[2]], b - 50)
    
    keepBurnin[b, ] <- params[-6]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[301:400, ]),
    stats::cov(keepBurnin[301:400, ])
  )
  
  # for 4
  for (b in 401:1000) {
    params <- rwBmetropolisSQAR1K1(params, sd * (Sigma[[2]] + I), T, L, n, Y, R.phi, Id, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-6], Sigma[[1]], Sigma[[2]], b - 300)
    
    keepBurnin[b, ] <- params[-6]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[751:1000, ]),
    stats::cov(keepBurnin[751:1000, ])
  )
  
  # for 5
  for (b in 1001:2000) {
    params <- rwBmetropolisSQAR1K1(params, sd * (Sigma[[2]] + I), T, L, n, Y, R.phi, Id, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-6], Sigma[[1]], Sigma[[2]], b - 750)
    
    keepBurnin[b, ] <- params[-6]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[1501:2000, ]),
    stats::cov(keepBurnin[1501:2000, ])
  )
  
  # for 6
  for (b in 2001:10000) {
    params <- rwBmetropolisSQAR1K1(params, sd * (Sigma[[2]] + I), T, L, n, Y, R.phi, Id, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-6], Sigma[[1]], Sigma[[2]], b - 1500)
    
    keepBurnin[b, ] <- params[-6]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[6001:10000, ]),
    stats::cov(keepBurnin[6001:10000, ])
  )
  
  # burn-in
  if (n.burnin < 1000) n.burnin <- 1000
  for (b in 10001:(10000 + n.burnin)) {
    params <- rwBmetropolisSQAR1K1(params, sd * (Sigma[[2]] + I), T, L, n, Y, R.phi, Id, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-6], Sigma[[1]], Sigma[[2]], b - 9950)
    
    keepBurnin[b, ] <- params[-6]
  }
  
  Sigma <- sd * (stats::cov(keepBurnin[10000 + (1:n.burnin), ]) + I)
  
  print("Burn-in finished")
  
  # FINAL
  for (b in 1:n.sims) {
    if (b %% n.report == 0) print(paste("Iteration:", b))
    
    params <- rwBmetropolisSQAR1K1(params, Sigma, T, L, n, Y, R.phi, Id, prior, tol)
    
    # Sigma <- MuSigmaUpdate(params[-6], Sigma[[1]], Sigma[[2]], b + 50)
    
    if (b %% n.thin == 0) {
      keep[b / n.thin, ] <- params[-6]
    }
  }
  
  keep <- exp(keep)
  keep[,5] <- keep[,5] / (1 + keep[,5])
  
  return(keep)
}

#' Spatial GP QAR(1) with K = 1
#' 
#' @importFrom stats cov
#' @importFrom stats dist
#' 
#' @description 
#' \deqn{\eta_1(s;\tau) = F_{a_1(s),b_1(s)}(\tau)}
#' \deqn{\eta_2(s;\tau) = F_{a_2(s),b_2(s)}(\tau)}
#' 
#' Burn-in of 100,000 iterations. Run at least 1,000,000 n.sims (might take 8/9 hours)
#' 
#' @param Y ARRAY, rows are independent, columns have autoregression, 
#'   third dim have copula dep (space)
#' @param coords COORDINATES \eqn{n \times 2} of the locations
#' @param n.sims,n.thin,n.burnin,n.report (i) Number of iterations after 
#'   burn-in. (ii) Thinning rate. (iii) Number of iterations discarded at the 
#'   beginning. (iv) Report the number of iterations rate.
#' @param inits Initial values
#' @param prior Vector (length 2). Standard deviation of the normal prior 
#'   for (1) \eqn{a}'s and \eqn{b}'s, and (2)
#'   \eqn{\log \sigma}'s, respectively
#' @param tol Tolerance in the univariate rootfinder. A value too small can lead 
#'   to numerical overflow in the likelihood 
#' @return MATRIX, rows are simulations, columns are parameters
#'   \deqn{\gamma,\exp a_1,\exp b_1,\exp a_2,\exp b_2,
#'   \sigma_{a_1}^2,\sigma_{b_1}^2,\sigma_{a_2}^2,\sigma_{b_2}^2,
#'   a_1(s_1),\ldots,a_1(s_n),b_1(s_1),\ldots,a_2(s_1),\ldots,b_2(s_n)}
#' @export 
SGPQAR1K1 <- function(Y, coords, 
                      n.sims = 1000000,
                      n.thin = 1000,
                      n.burnin = 100000,
                      n.report = 1000, 
                      inits = rep(0, 9 + 4 * dim(Y)[3]), 
                      prior = c(3, 3), 
                      tol = 10e-16) {
  
  T <- dim(Y)[1]
  L <- dim(Y)[2]
  n <- dim(Y)[3]
  d <- 9 + 4 * n
  d1 <- d + 1
  sd <- 2.4^2 / d
  epsilon <- 1e-06
  params <- inits
  keepBurnin <- matrix(nrow = 100000 + n.burnin, ncol = d)
  keep       <- matrix(nrow = n.sims / n.thin, ncol = d)
  colnames(keep) <- c("gamma", "a1", "b1", "a2", "b2", 
    "sigmaa1", "sigmab1", "sigmaa2", "sigmab2",
    do.call(paste0, expand.grid(1:n, c("a1s", "b1s", "a2s", "b2s"))[,2:1]))
  I <- epsilon * diag(d)
  
  if (2 * d + 1 > 1000) stop("Too many spatial locations (change code for pre burn-in)")
  
  distance <- stats::dist(coords)
  phi <- 3 / max(distance) # fixed phi
  dAux <- matrix(0, nrow = n, ncol = n)
  dAux[lower.tri(dAux)] <- distance
  dAux <- dAux + t(dAux)
  R.phi <- exp(- phi * dAux)
  Rinv <- solve(R.phi)
  Id <- diag(n)
  
  params <- c(params, logfallSGPQAR1K1(T, L, n, inits, Y, R.phi, Rinv, Id, prior, tol))
  
  # for 1
  for (b in 1:(2 * d)) {
    params <- rwBmetropolisSGPQAR1K1(params, sd * I, T, L, n, Y, R.phi, Rinv, Id, prior, tol)
    
    keepBurnin[b, ] <- params[-d1]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[1:(2 * d), ]),
    stats::cov(keepBurnin[1:(2 * d), ])
  )

  # for 2
  for (b in (2 * d + 1):1000) {
    params <- rwBmetropolisSGPQAR1K1(params, sd * (Sigma[[2]] + I), T, L, n, Y, R.phi, Rinv, Id, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-d1], Sigma[[1]], Sigma[[2]], b - 50)
    
    keepBurnin[b, ] <- params[-d1]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[501:1000, ]),
    stats::cov(keepBurnin[501:1000, ])
  )
  
  # for 3
  for (b in 1001:5000) {
    params <- rwBmetropolisSGPQAR1K1(params, sd * (Sigma[[2]] + I), T, L, n, Y, R.phi, Rinv, Id, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-d1], Sigma[[1]], Sigma[[2]], b - 500)
    
    keepBurnin[b, ] <- params[-d1]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[3001:5000, ]),
    stats::cov(keepBurnin[3001:5000, ])
  )
  
  # for 4
  for (b in 5001:10000) {
    params <- rwBmetropolisSGPQAR1K1(params, sd * (Sigma[[2]] + I), T, L, n, Y, R.phi, Rinv, Id, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-d1], Sigma[[1]], Sigma[[2]], b - 3000)
    
    keepBurnin[b, ] <- params[-d1]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[8001:10000, ]),
    stats::cov(keepBurnin[8001:10000, ])
  )
  
  # for 5
  for (b in 10001:50000) {
    params <- rwBmetropolisSGPQAR1K1(params, sd * (Sigma[[2]] + I), T, L, n, Y, R.phi, Rinv, Id, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-d1], Sigma[[1]], Sigma[[2]], b - 8000)
    
    keepBurnin[b, ] <- params[-d1]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[40001:50000, ]),
    stats::cov(keepBurnin[40001:50000, ])
  )
  
  # for 6
  for (b in 50001:100000) {
    params <- rwBmetropolisSGPQAR1K1(params, sd * (Sigma[[2]] + I), T, L, n, Y, R.phi, Rinv, Id, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-d1], Sigma[[1]], Sigma[[2]], b - 40000)
    
    keepBurnin[b, ] <- params[-d1]
  }
  
  Sigma <- list(
    colMeans(keepBurnin[90001:100000, ]),
    stats::cov(keepBurnin[90001:100000, ])
  )
  
  # burn-in
  if (n.burnin < 1000) n.burnin <- 1000
  for (b in 10001:(10000 + n.burnin)) {
    params <- rwBmetropolisSGPQAR1K1(params, sd * (Sigma[[2]] + I), T, L, n, Y, R.phi, Rinv, Id, prior, tol)
    
    Sigma <- MuSigmaUpdate(params[-d1], Sigma[[1]], Sigma[[2]], b - 99500)
    
    keepBurnin[b, ] <- params[-d1]
  }
  
  Sigma <- sd * (stats::cov(keepBurnin[100000 + (1:n.burnin), ]) + I)
  
  print("Burn-in finished")
  
  # FINAL
  for (b in 1:n.sims) {
    if (b %% n.report == 0) print(paste("Iteration:", b))
    
    params <- rwBmetropolisSGPQAR1K1(params, Sigma, T, L, n, Y, R.phi, Rinv, Id, prior, tol)
    
    # Sigma <- MuSigmaUpdate(params[-d1], Sigma[[1]], Sigma[[2]], b + 100)
    
    if (b %% n.thin == 0) {
      keep[b / n.thin, ] <- params[-d1]
    }
  }
  
  keep <- exp(keep)
  keep[,1] <- keep[,1] / (1 + keep[,1])
  
  return(keep)
}