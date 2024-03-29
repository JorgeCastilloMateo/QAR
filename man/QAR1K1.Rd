% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/QAR1.R
\name{QAR1K1}
\alias{QAR1K1}
\title{QAR(1) with K = 1}
\usage{
QAR1K1(
  Y,
  n.sims = 10000,
  n.thin = 10,
  n.burnin = 10000,
  n.report = 1000,
  inits = rep(0, 4),
  prior = 3,
  tol = 1e-15
)
}
\arguments{
\item{Y}{MATRIX, rows are independent, columns have autoregression}

\item{n.sims, n.thin, n.burnin, n.report}{(i) Number of iterations after 
burn-in. (ii) Thinning rate. (iii) Number of iterations discarded at the 
beginning. (iv) Report the number of iterations rate.}

\item{inits}{Initial values (transformed scale \eqn{(-\infty,\infty)})}

\item{prior}{Number. Standard deviation of the zero-mean normal prior for 
\eqn{\log a}'s and \eqn{\log b}'s}

\item{tol}{Tolerance in the univariate rootfinder. A value too small can lead 
to numerical overflow in the likelihood}
}
\value{
A \code{"QAR1K1"} list with elements:
  \item{params}{Matrix where rows are simulations and cols are parameters
  \deqn{a_1,b_1,a_2,b_2}}
  \item{\code{y}}{Data fitted}
}
\description{
This function fits the model 
\deqn{Q_{Y_{tl}}(\tau \mid y_{t,l-1}) = y_{t,l-1} \eta_1(\tau) + (1 - y_{t,l-1}) \eta_2(\tau)}
with
\deqn{\eta_1(\tau) = F(\tau \mid a_1,b_1)}
\deqn{\eta_2(\tau) = F(\tau \mid a_2,b_2)}
and \eqn{F(\tau \mid a,b)} the Kumaraswamy cdf with parameters \eqn{a} and \eqn{b}
}
