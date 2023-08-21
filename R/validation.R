#' Model Validation
#' 
#' @importFrom extraDistr pkumar
#' @importFrom stats qnorm
#' @importFrom stats quantile
#' 
#' @description Validate the predicted quantiles with
#' \deqn{p(\tau) = \frac{1}{T-1}\sum_{t=2}^{T} E[\textbf{1}(y_{t} < Q_{Y_{t}}(\tau \mid y_{t-1}))]}
#' and
#' \deqn{p = \sqrt{\int_{0}^{1} \left(\frac{p(\tau) - \tau}{\sqrt{(\tau - \tau^2) / (T-1)}}\right)^2 \,d\tau} \approx \sqrt{\frac{1}{\# G} \sum_{\tau \in G} \left(\frac{p(\tau) - \tau}{\sqrt{(\tau - \tau^2) / (T-1)}}\right)^2}}
#' with \eqn{G \subset (0,1)} a grid of values for \eqn{\tau}.
#' 
#' Also, 
#' \deqn{\delta(\tau) = \frac{1}{T-1}\sum_{t=2}^{T} \delta_{\tau}(y_{t} - Q_{Y_{t}}(\tau \mid y_{t-1}))}
#' where \eqn{\delta_{\tau}(u) = u(\tau - \textbf{1}(u<0))}, and
#' \deqn{\int_{0}^{1} \delta(\tau) d\tau \approx \frac{1}{\# G} \sum_{\tau \in G} \delta(\tau)}
#' with \eqn{G \subset (0,1)} as above. 
#' 
#' Also, \eqn{R^1(\tau)} and \eqn{R^1}.
#' 
#' @aliases validation.default validation.QAR1K1 validation.QAR1K2 
#'   validation.QAR1K validation.QAR2K1 validation.KX2006 validation
#'   
#' @param model Output from functions \code{\link{QAR1K1}}, 
#'   \code{\link{QAR1K2}}, \code{\link{QAR1K}}, \code{\link{QAR2K1}}, or
#'   \code{\link{KX2006}}
#'   (\code{y} in \code{model} must be VECTOR dimensional)
#' @param global.measure Logical. If \code{TRUE}, returns the measurements
#'   integrated, otherwise returns one measurement for each \eqn{\tau}
#' @param tau VECTOR of the grid of \eqn{\tau}-quantiles
#' @return Measurement of \eqn{p}, \eqn{\delta}, and \eqn{R^1} (VECTOR) or 
#'   measurements of \eqn{p(\tau)}, \eqn{\delta(\tau)}, and \eqn{R^1(\tau)}
#'   (MATRIX, rows \eqn{\tau}'s cols \eqn{p}, \eqn{\delta}, and \eqn{R^1}).
#' @export validation
validation <- function(model, global.measure = TRUE, tau = 1:99 / 100) {
  
  UseMethod("validation", model)
}

#' @rdname validation
#' @method validation default
#' @export validation.default
#' @export 
validation.default <- function(model, global.measure = TRUE, tau = 1:99 / 100) {
  stop("'class(model)' is not supported. Must be one of 'QAR1K1', 'QAR1K2', 'QAR1K', 'QAR2K1', 'KX2006'")
}

#' @rdname validation
#' @method validation QAR1K1
#' @export validation.QAR1K1
#' @export 
validation.QAR1K1 <- function(model, global.measure = TRUE, tau = 1:99 / 100) {
  
  y <- model$y
  T <- length(y)
  p <- matrix(nrow = length(tau), ncol = T - 1)
  d <- matrix(nrow = length(tau), ncol = T - 1)
  R <- matrix(nrow = length(tau), ncol = T - 1)
  Qmarg <- stats::quantile(y, tau)
  delta <- function(u, tau) {
    u * (tau - (u < 0))
  }
  for (q in 1:length(tau)) {
    eta1 <- extraDistr::pkumar(tau[q], model$params[,"a1"], model$params[,"b1"])
    eta2 <- extraDistr::pkumar(tau[q], model$params[,"a2"], model$params[,"b2"])
    for (t in 2:T) {
      Res <- y[t] - (y[t - 1] * eta1 + (1 - y[t - 1]) * eta2)
      p[q, t - 1] <- mean(Res < 0)
      d[q, t - 1] <- delta(mean(Res), tau[q])
      R[q, t - 1] <- delta(y[t] - Qmarg[q], tau[q])
    }
  }
  
  p <- rowMeans(p)
  d <- rowMeans(d)
  R <- 1 - d / rowMeans(R)
  
  if (global.measure) {
    p <- sqrt(mean((T - 1) * (p - tau)^2 / (tau - tau^2)))
    d <- mean(d)
    R <- mean(R)
    return(data.frame("p" = p, "delta" = d, "R" = R))
  } else {
    x <- data.frame("p" = p, "delta" = d, "R" = R)
    rownames(x) <- tau
    return(x)
  }
}

#' @rdname validation
#' @method validation QAR1K2
#' @export validation.QAR1K2
#' @export 
validation.QAR1K2 <- function(model, global.measure = TRUE, tau = 1:99 / 100) {
  
  y <- model$y
  T <- length(y)
  p <- matrix(nrow = length(tau), ncol = T - 1)
  d <- matrix(nrow = length(tau), ncol = T - 1)
  R <- matrix(nrow = length(tau), ncol = T - 1)
  Qmarg <- stats::quantile(y, tau)
  delta <- function(u, tau) {
    u * (tau - (u < 0))
  }
  for (q in 1:length(tau)) {
    eta1 <- model$params[,"lambda1"] * extraDistr::pkumar(tau[q], model$params[,"a1"], model$params[,"b1"]) + 
      (1 - model$params[,"lambda1"]) * extraDistr::pkumar(tau[q], model$params[,"a2"], model$params[,"b2"])
    eta2 <- model$params[,"lambda2"] * extraDistr::pkumar(tau[q], model$params[,"a3"], model$params[,"b3"]) + 
      (1 - model$params[,"lambda2"]) * extraDistr::pkumar(tau[q], model$params[,"a4"], model$params[,"b4"])
    for (t in 2:T) {
      Res <- y[t] - (y[t - 1] * eta1 + (1 - y[t - 1]) * eta2)
      p[q, t - 1] <- mean(Res < 0)
      d[q, t - 1] <- delta(mean(Res), tau[q])
      R[q, t - 1] <- delta(y[t] - Qmarg[q], tau[q])
    }
  }
  
  p <- rowMeans(p)
  d <- rowMeans(d)
  R <- 1 - d / rowMeans(R)
  
  if (global.measure) {
    p <- sqrt(mean((T - 1) * (p - tau)^2 / (tau - tau^2)))
    d <- mean(d)
    R <- mean(R)
    return(data.frame("p" = p, "delta" = d, "R" = R))
  } else {
    x <- data.frame("p" = p, "delta" = d, "R" = R)
    rownames(x) <- tau
    return(x)
  }
}

#' @rdname validation
#' @method validation QAR1K
#' @export validation.QAR1K
#' @export 
validation.QAR1K <- function(model, global.measure = TRUE, tau = 1:99 / 100) {
  
  y <- model$y
  coef <- model$coef
  
  K <- length(coef) / 2

  T <- length(y)
  p <- matrix(nrow = length(tau), ncol = T - 1)
  d <- matrix(nrow = length(tau), ncol = T - 1)
  R <- matrix(nrow = length(tau), ncol = T - 1)
  Qmarg <- stats::quantile(y, tau)
  delta <- function(u, tau) {
    u * (tau - (u < 0))
  }
  for (q in 1:length(tau)) {
    eta1 <- (1 - rowSums(model$params[,1:(K - 1)])) * 
      extraDistr::pkumar(tau[q], coef[2 * K - 1], coef[2 * K])
    eta2 <- (1 - rowSums(model$params[,K:(2 * K - 2)])) * 
      extraDistr::pkumar(tau[q], coef[2 * K - 1], coef[2 * K])
    for (k in 1:(K-1)) {
      eta1 <- eta1 + model$params[,k] * extraDistr::pkumar(tau[q], coef[2 * k - 1], coef[2 * k])
      eta2 <- eta2 + model$params[,k + K - 1] * extraDistr::pkumar(tau[q], coef[2 * k - 1], coef[2 * k])
    }
    
    for (t in 2:T) {
      Res <- y[t] - (y[t - 1] * eta1 + (1 - y[t - 1]) * eta2)
      p[q, t - 1] <- mean(Res < 0)
      d[q, t - 1] <- delta(mean(Res), tau[q])
      R[q, t - 1] <- delta(y[t] - Qmarg[q], tau[q])
    }
  }
  
  p <- rowMeans(p)
  d <- rowMeans(d)
  R <- 1 - d / rowMeans(R)
  
  if (global.measure) {
    p <- sqrt(mean((T - 1) * (p - tau)^2 / (tau - tau^2)))
    d <- mean(d)
    R <- mean(R)
    return(data.frame("p" = p, "delta" = d, "R" = R))
  } else {
    x <- data.frame("p" = p, "delta" = d, "R" = R)
    rownames(x) <- tau
    return(x)
  }
}

#' @rdname validation
#' @method validation QAR2K1
#' @export validation.QAR2K1
#' @export 
validation.QAR2K1 <- function(model, global.measure = TRUE, tau = 1:99 / 100) {
  
  y <- model$y
  T <- length(y)
  p <- matrix(nrow = length(tau), ncol = T - 2)
  d <- matrix(nrow = length(tau), ncol = T - 2)
  R <- matrix(nrow = length(tau), ncol = T - 2)
  Qmarg <- stats::quantile(y, tau)
  delta <- function(u, tau) {
    u * (tau - (u < 0))
  }
  for (q in 1:length(tau)) {
    eta1 <- extraDistr::pkumar(tau[q], model$params[,"a1"], model$params[,"b1"])
    eta2 <- extraDistr::pkumar(tau[q], model$params[,"a2"], model$params[,"b2"])
    eta3 <- extraDistr::pkumar(tau[q], model$params[,"a3"], model$params[,"b3"])
    for (t in 3:T) {
      Res <- y[t] - (model$params[,"pi"] * y[t - 1] * eta1 + (1 - model$params[,"pi"]) * y[t - 2] * eta2 + (1 - model$params[,"pi"] * y[t - 1] - (1 - model$params[,"pi"]) * y[t - 2]) * eta3)
      p[q, t - 2] <- mean(Res < 0)
      d[q, t - 2] <- delta(mean(Res), tau[q])
      R[q, t - 2] <- delta(y[t] - Qmarg[q], tau[q])
    }
  }
  
  p <- rowMeans(p)
  d <- rowMeans(d)
  R <- 1 - d / rowMeans(R)
  
  if (global.measure) {
    p <- sqrt(mean((T - 2) * (p - tau)^2 / (tau - tau^2)))
    d <- mean(d)
    R <- mean(R)
    return(data.frame("p" = p, "delta" = d, "R" = R))
  } else {
    x <- data.frame("p" = p, "delta" = d, "R" = R)
    rownames(x) <- tau
    return(x)
  }
}

#' @rdname validation
#' @method validation KX2006
#' @export 
validation.KX2006 <- function (model, global.measure = TRUE, tau = 1:99/100) {
  
  y <- model$y
  T <- length(y)
  p <- matrix(nrow = length(tau), ncol = T - 1)
  d <- matrix(nrow = length(tau), ncol = T - 1)
  R <- matrix(nrow = length(tau), ncol = T - 1)
  Qmarg <- stats::quantile(y, tau)
  delta <- function(u, tau) {
    u * (tau - (u < 0))
  }
  for (q in 1:length(tau)) {
    for (t in 2:T) {
      Res <- y[t] - (model$params[,"mu"] + model$params[,"sigma"] * stats::qnorm(tau[q]) + pmin(model$params[,"gamma0"] + model$params[,"gamma1"] * tau[q], 1) * y[t-1])
      p[q, t - 1] <- mean(Res < 0)
      d[q, t - 1] <- delta(mean(Res), tau[q])
      R[q, t - 1] <- delta(y[t] - Qmarg[q], tau[q])
    }
  }
  p <- rowMeans(p)
  d <- rowMeans(d)
  R <- 1 - d/rowMeans(R)
  if (global.measure) {
    p <- sqrt(mean((T - 1) * (p - tau)^2/(tau - tau^2)))
    d <- mean(d)
    R <- mean(R)
    return(data.frame("p" = p, "delta" = d, "R" = R))
  }
  else {
    x <- data.frame("p" = p, "delta" = d, "R" = R)
    rownames(x) <- tau
    return(x)
  }
}
