% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validation.R
\name{validation}
\alias{validation}
\alias{validation.default}
\alias{validation.QAR1K1}
\alias{validation.QAR1K2}
\alias{validation.QAR1K}
\alias{validation.QAR2K1}
\alias{validation.KX2006}
\title{Model Validation}
\usage{
validation(model, global.measure = TRUE, tau = 1:99/100)

\method{validation}{default}(model, global.measure = TRUE, tau = 1:99/100)

\method{validation}{QAR1K1}(model, global.measure = TRUE, tau = 1:99/100)

\method{validation}{QAR1K2}(model, global.measure = TRUE, tau = 1:99/100)

\method{validation}{QAR1K}(model, global.measure = TRUE, tau = 1:99/100)

\method{validation}{QAR2K1}(model, global.measure = TRUE, tau = 1:99/100)

\method{validation}{KX2006}(model, global.measure = TRUE, tau = 1:99/100)
}
\arguments{
\item{model}{Output from functions \code{\link{QAR1K1}}, 
\code{\link{QAR1K2}}, \code{\link{QAR1K}}, \code{\link{QAR2K1}}, or
\code{\link{KX2006}}
(\code{y} in \code{model} must be VECTOR dimensional)}

\item{global.measure}{Logical. If \code{TRUE}, returns the measurements
integrated, otherwise returns one measurement for each \eqn{\tau}}

\item{tau}{VECTOR of the grid of \eqn{\tau}-quantiles}
}
\value{
Measurement of \eqn{p}, \eqn{\delta}, and \eqn{R^1} (VECTOR) or 
  measurements of \eqn{p(\tau)}, \eqn{\delta(\tau)}, and \eqn{R^1(\tau)}
  (MATRIX, rows \eqn{\tau}'s cols \eqn{p}, \eqn{\delta}, and \eqn{R^1}).
}
\description{
Validate the predicted quantiles with
\deqn{p(\tau) = \frac{1}{T-1}\sum_{t=2}^{T} E[\textbf{1}(y_{t} < Q_{Y_{t}}(\tau \mid y_{t-1}))]}
and
\deqn{p = \sqrt{\int_{0}^{1} \left(\frac{p(\tau) - \tau}{\sqrt{(\tau - \tau^2) / (T-1)}}\right)^2 \,d\tau} \approx \sqrt{\frac{1}{\# G} \sum_{\tau \in G} \left(\frac{p(\tau) - \tau}{\sqrt{(\tau - \tau^2) / (T-1)}}\right)^2}}
with \eqn{G \subset (0,1)} a grid of values for \eqn{\tau}.

Also, 
\deqn{\delta(\tau) = \frac{1}{T-1}\sum_{t=2}^{T} \delta_{\tau}(y_{t} - Q_{Y_{t}}(\tau \mid y_{t-1}))}
where \eqn{\delta_{\tau}(u) = u(\tau - \textbf{1}(u<0))}, and
\deqn{\int_{0}^{1} \delta(\tau) d\tau \approx \frac{1}{\# G} \sum_{\tau \in G} \delta(\tau)}
with \eqn{G \subset (0,1)} as above. 

Also, \eqn{R^1(\tau)} and \eqn{R^1}.
}
