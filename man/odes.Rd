% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/odes.R
\name{odes}
\alias{odes}
\title{Observation-driven Exponential Smoothing}
\usage{
odes(
  formula,
  data,
  subset,
  offset,
  contrasts = NULL,
  Huber = TRUE,
  kappa = 1.345,
  max_lag = NULL,
  ...
)
}
\arguments{
\item{formula}{Symbol description of the regression model of type \code{y ~ x}.}

\item{data}{An optional data frame, list or environment containing the
variables in the regression model.}

\item{subset}{An optional vector specifying a subset of observations to
be used in the fitting process.}

\item{offset}{One or more \link{offset} terms can be included in the formula
instead or as well, and if more than one is specified their sum is used.
See \link{model.offset}.}

\item{contrasts}{An optional list. See the contrasts.arg of
\link{model.matrix.default}.}

\item{Huber}{Logical. Should "Huber-psi" robustification be applied before
exponential smoothing?}

\item{kappa}{Numeric. The tuning parameter of the Huber psi function that is
used if \code{Huber = TRUE}.}

\item{max_lag}{Numeric. Maximal number of lags used in the HAC sandwich
variance. If \code{NULL} then is computed as \code{floor(0.75 * n_obs ^ (0.33))}.}

\item{\\dots}{Additional optional parameters, currently only the \code{control}
argument of \link{optim}.}
}
\value{
An object of class \code{odes}.
}
\description{
Computes exponential smoothing for Poisson counts using the approach
described in Karlis, Pedeli and Varin (2023). \code{odes} assumes that the
distribution of \eqn{Y_t} given the past history is
\eqn{Y_t| \mathcal{F}_{t-1} \sim \text{Poisson}(\mu_t)} with the state
variable \eqn{W_t=\log \mu_t} defined recursively as
\deqn{
\begin{aligned}
W_t &= x_t^T \beta + Z_t, \\
Z_t &= \alpha \epsilon_{t-1} + (1-\alpha) Z_{t-1},
\end{aligned}
}
where \eqn{x_t} is a vector of covariates, \eqn{\beta} is the corresponding
vector of regression coefficients, \eqn{\alpha \in [0,1]} is a smoothing
parameter and \eqn{\epsilon_t=(Y_t-e^{w_t})/e^{w_t/2}} are the Pearson
residuals. This form of exponential smoothing is a special case of the
general class of observation-driven models described in Davis, Dunsmuir
and Streett (2003).
}
\examples{
dengue_data <- data.frame(counts = dengue, time = seq_len(length(dengue)))
model <- counts ~ sin(2 * pi * time / 12) + cos(2 * pi * time / 12)
## non-robust odes
fit1 <- odes(model, data = dengue_data, Huber = FALSE)
fit1
## robust odes
fit2 <- odes(model, data = dengue_data)
fit2

}
\references{
{
Davis, R. A., Dunsmuir, W. T. and Streett, S. B. (2003). Observation‐driven
models for Poisson counts. \emph{Biometrika}, \strong{90} (4), 777-790.

Karlis, D., Pedeli, X. and Varin, C. (2023). Observation-driven exponential
smoothing for Poisson counts. \emph{Stat}, \strong{12} (1), e642.
}
}
