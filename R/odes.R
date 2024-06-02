#' Observation-driven Exponential Smoothing
#'
#' Computes exponential smoothing for Poisson counts using the approach
#' described in Karlis, Pedeli and Varin (2023). `odes` assumes that the
#' distribution of \eqn{Y_t} given the past history is
#' \eqn{Y_t| \mathcal{F}_{t-1} \sim \text{Poisson}(\mu_t)} with the state
#' variable \eqn{W_t=\log \mu_t} defined recursively as
#' \deqn{
#' \begin{aligned}
#' W_t &= x_t^T \beta + Z_t, \\
#' Z_t &= \alpha \epsilon_{t-1} + (1-\alpha) Z_{t-1},
#' \end{aligned}
#' }
#' where \eqn{x_t} is a vector of covariates, \eqn{\beta} is the corresponding
#' vector of regression coefficients, \eqn{\alpha \in [0,1]} is a smoothing
#' parameter and \eqn{\epsilon_t=(Y_t-e^{w_t})/e^{w_t/2}} are the Pearson
#' residuals. This form of exponential smoothing is a special case of the
#' general class of observation-driven models described in Davis, Dunsmuir
#' and Streett (2003).
#'
#' @param formula Symbol description of the regression model of type `y ~ x`.
#' @param data An optional data frame, list or environment containing the
#' variables in the regression model.
#' @param subset An optional vector specifying a subset of observations to
#' be used in the fitting process.
#' @param offset One or more [offset] terms can be included in the formula
#' instead or as well, and if more than one is specified their sum is used.
#' See [model.offset].
#' @param contrasts An optional list. See the contrasts.arg of
#' [model.matrix.default].
#' @param Huber Logical. Should "Huber-psi" robustification be applied before
#' exponential smoothing?
#' @param kappa Numeric. The tuning parameter of the Huber psi function that is
#' used if `Huber = TRUE`.
#' @param max_lag Numeric. Maximal number of lags used in the HAC sandwich
#' variance. If `NULL` then is computed as `floor(0.75 * n_obs ^ (0.33))`.
#' @param \\dots Additional optional parameters, currently only the `control`
#' argument of [optim].
#'
#' @return An object of class `odes`.
#'
#' @examples
#' dengue_data <- data.frame(counts = dengue, time = seq_len(length(dengue)))
#' model <- counts ~ sin(2 * pi * time / 12) + cos(2 * pi * time / 12)
#' ## non-robust odes
#' fit1 <- odes(model, data = dengue_data, Huber = FALSE)
#' fit1
#' ## robust odes
#' fit2 <- odes(model, data = dengue_data)
#' fit2
#'
#' @references{
#' Davis, R. A., Dunsmuir, W. T. and Streett, S. B. (2003). Observation‐driven
#' models for Poisson counts. \emph{Biometrika}, \strong{90} (4), 777-790.
#'
#' Karlis, D., Pedeli, X. and Varin, C. (2023). Observation-driven exponential
#' smoothing for Poisson counts. \emph{Stat}, \strong{12} (1), e642.
#' }
#'
#' @export

odes <- function(formula, data, subset, offset, contrasts = NULL, Huber = TRUE,
                 kappa = 1.345, max_lag = NULL, ...)
{

  ## lines below are inherited from glm
  call <- match.call()
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = TRUE)
  m <- match(c("formula", "data", "subset", "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  y <- model.response(mf, "any")
  x <- if (!is.empty.model(mt))
    model.matrix(mt, mf, contrasts)
  else matrix(NA, length(y), 0L)
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != length(y))
      stop(gettextf("number of offsets is %d should equal %d (number of
                    observations)", length(offset), NROW(y)), domain = NA)
  }
  if (is.null(offset))
    offset <- rep.int(0, length(y))
  ## ends lines inherited from glm
  if (!is.numeric(kappa) || kappa <= 0)
    stop("'kappa' must be positive")
  n_obs <- length(y)
  p <- NCOL(x)
  if (is.null(max_lag))
    max_lag <- floor(0.75 * n_obs ^ (0.33))
  if (!is.numeric(max_lag) || max_lag <= 0)
    stop("'max_lag' must be positive")
  ## create list with results
  ans <- structure(list(n_obs = n_obs, p = p, y = y, x = x, offset = offset,
                        data = data, call = match.call(), terms = mt,
                        Huber = Huber, kappa = kappa, max_lag = max_lag),
                   class = "odes")

  ## 1. Starting values
  if (ans$Huber) {
    mod0 <- robustbase::glmrob(eval(ans$call$formula), data = data,
                               family = poisson)
    if (is.null(ans$kappa))
      ans$kappa <- mod0$tcc
    res <- (y - fitted(mod0)) / sqrt(fitted(mod0))
    bad_obs <- (abs(res) >= ans$kappa)
    res[bad_obs] <- ans$kappa * sign(res[bad_obs])
    ## adjusted responses
    adj_y <- res * sqrt(fitted(mod0)) + fitted(mod0)
  }
  else {
    mod0 <- glm.fit(x, y, offset = offset, family = poisson())
    adj_y <- y
  }
  ans$adj_y <- adj_y
  ans$beta_glm <- coef(mod0)

  ## 2. Fitting
  ## negative log likelihood
  ans$neglik <- function(beta, alpha) {
    eta <- x %*% beta + offset
    compute_loglik(eta, adj_y, alpha)$nlik
  }
  ## vector of individual scores
  ans$scores <- function(beta, alpha) {
    eta <- x %*% beta + offset
    compute_scores(eta, x, adj_y, alpha)
  }
  lik_score <- function(beta, alpha) colSums(ans$scores(beta, alpha))
  ## ... and go!
  alphas <- seq(0, .99, by = .01)
  objfun_alpha <- rep(NA, length(alphas))
  betas <- matrix(nrow = p, ncol = length(alphas))
  beta_start <- ans$beta_glm
  rownames(betas) <- names(beta_start)
  for (i in seq_len(length(alphas))) {
    op <- try(optim(beta_start, fn = ans$neglik, gr = lik_score,
                    alpha = alphas[i], method = "BFGS", ...), silent = TRUE)
    if (!inherits(op, "try-error")) {
      ## update
      betas[, i] <- op$par
      beta_start <- op$par
      objfun_alpha[i] <- op$value
    }
  }
  id_best_alpha <- which.min(objfun_alpha)
  ans$objfun_alpha <- objfun_alpha
  ans$alphas <- alphas
  ans$beta <- betas[, id_best_alpha]
  ans$alpha <- alphas[id_best_alpha]
  ans$param <- c(ans$beta, ans$alpha)
  ## extract model components
  eta <- x %*% ans$beta + offset
  fit <- compute_loglik(eta, adj_y, ans$alpha)
  ans$z <- fit$z
  ans$eps <- fit$eps
  ans$w <- fit$w
  ## compute vcov
  mean_grad <- function(x) colMeans(ans$scores(x, ans$alpha))
  ans$J <- try(numDeriv::jacobian(mean_grad, ans$beta), silent = TRUE)
  if (inherits(ans$J, "try-error"))
    stop("\nCannot compute the Hessian.")
  Jinv <- try(solve(ans$J), silent = TRUE)
  if (inherits(Jinv, "try-error"))
    stop("\nCannot invert the Hessian.")
  ## the sandwich filling
  ans$efun <- ans$scores(ans$beta, ans$alpha)
  acf.es <- acf(ans$efun, lag.max = max_lag, type = "cov", plot = FALSE,
                demean = FALSE)[[1]]
  vmat <- acf.es[1L, , ]
  autocov <- acf.es[-1L, , ]
  wb <- (1 - (1:max_lag) / max_lag)
  if (is.array(autocov))
    wsum <- apply(autocov, c(2, 3), function(x) sum(x * wb))
  else
    wsum <- sum(autocov * wb)
  K <- vmat + wsum + t(wsum)
  ans$vcov <- matrix(NA, nrow = p, ncol = p)
  ans$vcov <- Jinv / n_obs
  rownames(ans$vcov) <- colnames(ans$vcov) <- names(ans$beta)
  ans$vcov_robust <- matrix(NA, nrow = p, ncol = p)
  ans$vcov_robust <- try(Jinv %*% K %*% Jinv / n_obs, silent = TRUE)
  if (inherits(ans$vcov_robust, "try-error"))
    stop("\nCannot compute the sandwich variance.")
  rownames(ans$vcov_robust) <- colnames(ans$vcov_robust) <- names(ans$beta)
  ## bye bye
  class(ans) <- "odes"
  ans
}

#' Summarizing Observation-driven Exponential Smoothing
#'
#' `print` method for `odes` objects.
#'
#' @param x An object of `odes` class.
#' @param digits The number of significant digits to use when printing.
#' @param \\dots Not used.
#' @examples
#' dengue_data <- data.frame(counts = dengue, time = seq_len(length(dengue)))
#' fit <- odes(counts ~ sin(2 * pi * time / 12) + cos(2 * pi * time / 12),
#' data = dengue_data)
#' fit
#'
#' @references{
#' Karlis, D., Pedeli, X. and Varin, C. (2023). Observation-driven exponential
#' smoothing for Poisson counts. \emph{Stat}, \strong{12} (1), e642.
#' }
#'
#' @export
print.odes <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:", deparse(x$call, width.cutoff =
                           floor(getOption("width") * 0.85)), "", sep = "\n")
  cat("Regression parameters:\n")
  print.default(format(x$beta, digits = digits), print.gap = 2, quote = FALSE)
  cat("\n")
  cat("Smoothing parameter (alpha): ")
  cat(format(x$alpha, digits = digits))
  cat("\n")
  invisible(x)
}

#' Summarizing observation driven exponential smoothing
#'
#' `summary` method for `odes` objects.
#'
#' @param object An object of `odes` class.
#' @param sandwich Logical. Should standard errors be computed with the robust
#' sandwich HAC variance or not?
#' @param digits The number of significant digits to use when printing.
#' @param signif.stars Logical. If 'TRUE', `significance stars` are printed for
#' each coefficient.
#' @param \\dots Not used.
#'
#'#' @references{
#' Karlis, D., Pedeli, X. and Varin, C. (2023). Observation-driven exponential
#' smoothing for Poisson counts. \emph{Stat}, \strong{12} (1), e642.
#' }
#'
#' @examples
#' dengue_data <- data.frame(counts = dengue, time = seq_len(length(dengue)))
#' fit <- odes(counts ~ sin(2 * pi * time / 12) + cos(2 * pi * time / 12),
#' data = dengue_data)
#' ## summary with ordinary standard errors
#' summary(fit, sandwich = FALSE)
#' ## summary with robust sandwich standard errors
#' summary(fit)
#'
#' @export
summary.odes <- function (object, sandwich = TRUE,
                          digits = max(3L, getOption("digits") - 3L),
                          signif.stars = getOption("show.signif.stars"), ...)
{
  cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  beta_table <- object$beta
  if (sandwich)
    se_beta <- sqrt(diag(object$vcov_robust))
  else
    se_beta <- sqrt(diag(object$vcov))
  beta_table <- cbind(beta_table, se_beta, beta_table / se_beta,
                      2 * pnorm(-abs(beta_table / se_beta)))
  colnames(beta_table) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  rownames(beta_table) <- names(object$beta)
  printCoefmat(beta_table, digits = digits, signif.stars = signif.stars,
               na.print = "NA", ...)
  cat("\nSmoothing parameter (alpha):", object$alpha)
}

#' Prediction with Observation-driven Exponential Smoothing
#'
#' Computes one-step ahead forecasts with `odes`.
#'
#' @references{
#' Karlis, D., Pedeli, X. and Varin, C. (2023). Observation-driven exponential
#' smoothing for Poisson counts. \emph{Stat}, \strong{12} (1), e642.
#' }
#'
#' @param object An object of `odes` class.
#' @param newdata A data frame in which to look for variables with which to
#' predict.
#'
#' @examples
#' dengue_data <- data.frame(counts = dengue, time = seq_len(length(dengue)))
#' fit <- odes(counts ~ sin(2 * pi * time / 12) + cos(2 * pi * time / 12),
#' data = dengue_data[1:108, ])
#' predict(fit, dengue_data[109, ])
#' ## the true value is
#' dengue_data[109, "counts"]
#'
#' @export
predict.odes <- function(object, newdata) {
  if (!inherits(object, "odes"))
    stop("\nThe class of object must be 'odes'.\n")
  if (NROW(newdata) != 1)
    stop("\nBad specification of 'newdata'.\n")
  last_z <- object$z[object$n_obs]
  last_eps <- object$eps[object$n_obs]
  new_z <-  (1 - object$alpha) * last_z + object$alpha * last_eps
  ## next line needs to be tested with contrasts
  cleaned_terms <- delete.response(terms(object))
  new_x <- model.matrix(cleaned_terms, newdata)
  new_offset <- model.offset(model.frame(cleaned_terms, newdata))
  if (is.null(new_offset))
    new_offset <- 0
  ans <- exp(sum(new_x * object$beta) + new_offset + new_z)
  unname(ans)
}

#' Surveillance with Observation-driven Exponential Smoothing
#'
#' Monitors a time series of counts with observation-driven exponential smoothing.
#' Monitoring is carried out by comparing new observations (phase II data) with
#' an upper control limit calculated with parameters estimated from historical
#' data (phase I).
#'
#' @references{
#' Karlis, D., Pedeli, X. and Varin, C. (2023). Observation-driven exponential
#' smoothing for Poisson counts. \emph{Stat}, \strong{12} (1), e642.
#' }
#'
#' @param object An object of `odes` class corresponding to a fitted model to
#' phase I data.
#' @param new_data Data frame containing the observations to monitor
#' corresponding to phase II data.
#' @param level Level used in the computation the upper control limit. Default
#' is `0.99`.
#' @param plot Logical. Should a graph with the series and the upper control
#' limit be plotted? Default is `TRUE`.
#' @param names_as_dates Logical. Are the names of the observations dates?
#' Default is `FALSE`. Only used if `plot` is `TRUE`.
#' @param x_label Character. Label for x-axis. Default is `Time`. Only used if
#' `plot` is `TRUE`.
#' @param y_label Character. Label for y-axis. Default is  `Count`. Only used if
#' `plot` is `TRUE`.
#' @param \\dots Additional optional graphical parameters passed to `ggplot`.
#'
#' @return An object of class `surveillance.odes`.
#'
#' @examples
#' ## Replicates Section 4 of Karlis, Pedeli and Varin (2023)
#' dengue_data <- data.frame(counts = dengue, time = seq_len(length(dengue)))
#' dates <- seq(as.Date("2008/1/1"), by = "month", length.out =
#' nrow(dengue_data))
#' rownames(dengue_data) <- dates
#' ## phase I
#' dengue_phaseI <- dengue_data[1:108, ]
#' fit <- odes(counts ~ sin(2 * pi * time / 12) + cos(2 * pi * time / 12),
#' data = dengue_phaseI)
#' summary(fit)
#' ## phase 2
#' dengue_phaseII <- dengue_data[-(1:108), ]
#' surv <- surveillance.odes(fit, new_data = dengue_phaseII, names_as_dates =
#' TRUE, date_labels = "%Y", breaks = "2 year", x_label = NULL, y_label =
#' "Dengue cases")
#' surv$alarms
#'
#' @export
surveillance.odes <- function(object, new_data, level = 0.99, plot = TRUE,
                              names_as_dates = FALSE, x_label = "Time",
                              y_label = "Count", ...)
{

  if (!inherits(object, "odes"))
    stop("\nThe class of object must be 'odes'.\n")
  new_mf <- model.frame(object$terms, new_data)
  new_y <- model.response(new_mf)
  new_x <- model.matrix(new_mf, new_data)
  new_offset <- model.offset(new_mf)
  if (is.null(new_offset))
   new_offset <- rep.int(0, length(new_y))
  new_eta <- new_x %*% object$beta + new_offset
  if (object$Huber) {
    new_res <- (new_y - exp(new_eta)) / sqrt(exp(new_eta))
    new_bad_obs <- (abs(new_res) >= object$kappa)
    new_res[new_bad_obs] <- object$kappa * sign(new_res[new_bad_obs])
    ## adjusted responses
    new_adj_y <- new_res * sqrt(exp(new_eta)) + exp(new_eta)
    new_w <- compute_loglik(eta = new_eta, y = new_adj_y, alpha = object$alpha)$w
  }
  else
    new_w <- compute_loglik(eta = new_eta, y = new_y, alpha = object$alpha)$w
  new_upper <- qpois(level, exp(new_w))
  ans <- list(new_y = new_y, new_upper = new_upper,
              alarms = which(new_y > new_upper, arr.ind = TRUE))
  if (plot) {
    upper <- qpois(level, exp(object$w))
    up <- c(upper, new_upper)
    new_up <- c(rep(NA, length(object$y)), new_upper)
    if (names_as_dates)
      times <- as.Date(names(c(object$y, new_y)))
    else
      times <- seq_len(length(object$y) + length(new_y))
    gg_df <- data.frame(y = c(object$y, new_y), times = times)
    ans$gg <- ggplot(data = gg_df, aes(gg_df$times, gg_df$y)) +
      geom_line(color = gray(.5), size = 0.7) +
      geom_ribbon(aes(ymin = 0, ymax = up), fill = gray(.5), alpha = 0.3) +
      geom_ribbon(aes(ymin = 0, ymax = new_up), fill = "#D55E00D0",
                  alpha = 0.3) + ylab(y_label) + xlab(x_label) +
      theme_minimal() +
      theme(axis.title.y = element_text(angle = 0))

    if (names_as_dates)
      ans$gg <- ans$gg + scale_x_date(...)

    print(ans$gg)
  }
  ## bye bye
  class(ans) <- "surveillance.odes"
  invisible(ans)
}
