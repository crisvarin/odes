% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/odes.R
\name{predict.odes}
\alias{predict.odes}
\title{Prediction with Observation-driven Exponential Smoothing}
\usage{
\method{predict}{odes}(object, newdata)
}
\arguments{
\item{object}{An object of \code{odes} class.}

\item{newdata}{A data frame in which to look for variables with which to
predict.}
}
\description{
Computes one-step ahead forecasts with \code{odes}.
}
\examples{
dengue_data <- data.frame(counts = dengue, time = seq_len(length(dengue)))
fit <- odes(counts ~ sin(2 * pi * time / 12) + cos(2 * pi * time / 12),
data = dengue_data[1:108, ])
predict(fit, dengue_data[109, ])
## the true value is
dengue_data[109, "counts"]

}
\references{
{
Karlis, D., Pedeli, X. and Varin, C. (2023). Observation-driven exponential
smoothing for Poisson counts. \emph{Stat}, \strong{12} (1), e642.
}
}
