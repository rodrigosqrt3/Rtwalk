#' @importFrom stats sd quantile acf density fft nextn qnorm median var
#' @importFrom graphics par layout hist lines polygon
#' @importFrom grDevices adjustcolor
NULL

# Convergence diagnostics of Vehtari, Gelman, Simpson, Carpenter & Buerkner
# (2021), "Rank-normalization, folding, and localization: An improved R-hat
# for assessing convergence of MCMC", Bayesian Analysis 16(2), 667-718.
# The estimators reproduce those of the 'posterior' package; the test suite
# compares both. Throughout, `draws` is an iterations x chains matrix.

.is_constant_draws <- function(draws) {
  anyNA(draws) || any(!is.finite(draws)) ||
    abs(max(draws) - min(draws)) < .Machine$double.eps
}

.split_chains <- function(draws) {
  draws <- as.matrix(draws)
  n <- nrow(draws)
  half <- n / 2
  cbind(
    draws[seq_len(floor(half)), , drop = FALSE],
    draws[seq.int(ceiling(half + 1), n), , drop = FALSE]
  )
}

.z_scale <- function(draws) {
  ranks <- rank(draws, ties.method = "average")
  c <- 3 / 8
  z <- stats::qnorm((ranks - c) / (length(ranks) - 2 * c + 1))
  array(z, dim = dim(draws))
}

# Biased autocovariance estimate (Geyer, 1992) computed with the FFT.
.autocovariance <- function(x) {
  n <- length(x)
  var_x <- stats::var(x)
  if (var_x == 0) {
    return(rep(0, n))
  }
  m2 <- 2 * stats::nextn(n)
  centred <- c(x - mean(x), rep.int(0, m2 - n))
  ac <- Re(stats::fft(abs(stats::fft(centred))^2, inverse = TRUE)[seq_len(n)])
  ac / ac[1L] * var_x * (n - 1) / n
}

# Multi-chain effective sample size with Geyer's initial monotone sequence.
.ess_draws <- function(draws) {
  draws <- as.matrix(draws)
  n_chains <- ncol(draws)
  n <- nrow(draws)
  if (n < 3L || .is_constant_draws(draws)) {
    return(NA_real_)
  }

  acov <- apply(draws, 2, .autocovariance)
  acov_means <- rowMeans(acov)
  mean_var <- acov_means[1L] * n / (n - 1)
  var_plus <- mean_var * (n - 1) / n
  if (n_chains > 1L) {
    var_plus <- var_plus + stats::var(colMeans(draws))
  }

  rho_hat <- rep.int(0, n)
  t <- 0
  rho_hat_even <- 1
  rho_hat[t + 1] <- rho_hat_even
  rho_hat_odd <- 1 - (mean_var - acov_means[t + 2]) / var_plus
  rho_hat[t + 2] <- rho_hat_odd

  while (t < n - 5 && !is.nan(rho_hat_even + rho_hat_odd) &&
         (rho_hat_even + rho_hat_odd > 0)) {
    t <- t + 2
    rho_hat_even <- 1 - (mean_var - acov_means[t + 1]) / var_plus
    rho_hat_odd <- 1 - (mean_var - acov_means[t + 2]) / var_plus
    if ((rho_hat_even + rho_hat_odd) >= 0) {
      rho_hat[t + 1] <- rho_hat_even
      rho_hat[t + 2] <- rho_hat_odd
    }
  }
  max_t <- t
  if (rho_hat_even > 0) {
    rho_hat[max_t + 1] <- rho_hat_even
  }

  # Geyer's initial monotone sequence
  t <- 0
  while (t <= max_t - 4) {
    t <- t + 2
    if (rho_hat[t + 1] + rho_hat[t + 2] > rho_hat[t - 1] + rho_hat[t]) {
      rho_hat[t + 1] <- (rho_hat[t - 1] + rho_hat[t]) / 2
      rho_hat[t + 2] <- rho_hat[t + 1]
    }
  }

  ess <- n_chains * n
  tau_hat <- -1 + 2 * sum(rho_hat[seq_len(max_t)]) + rho_hat[max_t + 1]
  tau_hat <- max(tau_hat, 1 / log10(ess))
  ess / tau_hat
}

.rhat_draws <- function(draws) {
  draws <- as.matrix(draws)
  if (.is_constant_draws(draws)) {
    return(NA_real_)
  }
  n <- nrow(draws)
  chain_means <- colMeans(draws)
  chain_vars <- apply(draws, 2, stats::var)
  var_between <- n * stats::var(chain_means)
  var_within <- mean(chain_vars)
  sqrt((var_between / var_within + n - 1) / n)
}

# Rank-normalized split-R-hat: the maximum of the bulk and tail (folded)
# versions. Defined for a single chain as well, through splitting.
.rhat_rank <- function(draws) {
  if (.is_constant_draws(draws)) {
    return(NA_real_)
  }
  folded <- abs(draws - stats::median(draws))
  max(
    .rhat_draws(.z_scale(.split_chains(draws))),
    .rhat_draws(.z_scale(.split_chains(folded)))
  )
}

.ess_basic <- function(draws) .ess_draws(.split_chains(draws))

.ess_bulk <- function(draws) {
  if (.is_constant_draws(draws)) {
    return(NA_real_)
  }
  .ess_draws(.z_scale(.split_chains(draws)))
}

.ess_tail <- function(draws) {
  if (.is_constant_draws(draws)) {
    return(NA_real_)
  }
  ess_quantile <- function(prob) {
    indicator <- draws <= stats::quantile(draws, prob)
    storage.mode(indicator) <- "double"
    .ess_draws(.split_chains(indicator))
  }
  min(ess_quantile(0.05), ess_quantile(0.95))
}

# iterations x chains matrix of one parameter, truncating chains to a common
# length.
.parameter_draws <- function(chains, parameter_index) {
  common_length <- min(vapply(chains, nrow, integer(1)))
  vapply(
    chains,
    function(chain) chain[seq_len(common_length), parameter_index],
    numeric(common_length)
  )
}

.calculate_rhat <- function(chains) {
  vapply(
    seq_len(ncol(chains[[1L]])),
    function(j) {
      draws <- .parameter_draws(chains, j)
      if (nrow(as.matrix(draws)) < 4L) NA_real_ else .rhat_rank(draws)
    },
    numeric(1)
  )
}

.prepare_mcmc_chains <- function(samples, burnin_frac) {
  if (!is.numeric(burnin_frac) || length(burnin_frac) != 1L ||
      is.na(burnin_frac) || !is.finite(burnin_frac) ||
      burnin_frac < 0 || burnin_frac >= 1) {
    stop("`burnin_frac` must be one finite number in [0, 1).", call. = FALSE)
  }

  sample_chains <- if (is.list(samples) && !is.data.frame(samples)) {
    samples
  } else {
    list(samples)
  }

  if (length(sample_chains) == 0L) {
    stop("`samples` must contain at least one chain.", call. = FALSE)
  }

  sample_chains <- lapply(sample_chains, function(chain) {
    if (!is.matrix(chain) && !is.data.frame(chain)) {
      stop("Each chain in `samples` must be a numeric matrix or data frame.", call. = FALSE)
    }

    chain <- as.matrix(chain)
    if (!is.numeric(chain) || nrow(chain) < 2L || ncol(chain) < 1L ||
        any(!is.finite(chain))) {
      stop("Each chain must have at least two rows, one column, and only finite numeric values.", call. = FALSE)
    }

    chain
  })

  n_columns <- vapply(sample_chains, ncol, integer(1))
  if (length(unique(n_columns)) != 1L) {
    stop("All chains in `samples` must have the same number of columns.", call. = FALSE)
  }

  lapply(sample_chains, function(chain) {
    n_iter <- nrow(chain)
    n_burnin <- floor(burnin_frac * n_iter)
    keep <- seq.int(n_burnin + 1L, n_iter)

    if (length(keep) < 2L) {
      stop("`burnin_frac` must leave at least two samples in every chain.", call. = FALSE)
    }

    chain[keep, , drop = FALSE]
  })
}

.diagnostics_table <- function(post_burnin_chains, param_names = NULL) {
  post_burnin_samples <- do.call(rbind, post_burnin_chains)
  n_param <- ncol(post_burnin_samples)

  if (is.null(param_names)) {
    param_names <- colnames(post_burnin_samples)
    if (is.null(param_names) || anyNA(param_names) || any(param_names == "")) {
      param_names <- paste0("theta", seq_len(n_param))
    }
  } else if (!is.character(param_names) || length(param_names) != n_param ||
             anyNA(param_names)) {
    stop("`param_names` must be a character vector with one name per parameter.", call. = FALSE)
  }

  means <- colMeans(post_burnin_samples)
  sds <- apply(post_burnin_samples, 2, sd)
  quantiles <- apply(post_burnin_samples, 2, quantile, probs = c(0.025, 0.5, 0.975))

  per_parameter <- vapply(seq_len(n_param), function(j) {
    draws <- .parameter_draws(post_burnin_chains, j)
    if (nrow(as.matrix(draws)) < 4L) {
      return(rep(NA_real_, 5L))
    }
    ess <- .ess_basic(draws)
    c(
      ess = ess,
      ess_bulk = .ess_bulk(draws),
      ess_tail = .ess_tail(draws),
      iat = length(draws) / ess,
      rhat = .rhat_rank(draws)
    )
  }, numeric(5))

  ess <- per_parameter[1L, ]
  mcse <- ifelse(is.finite(ess) & ess > 0, sds / sqrt(ess), NA_real_)

  data.frame(
    Parameter = param_names,
    Mean = round(means, 4),
    SD = round(sds, 4),
    Q2.5 = round(quantiles["2.5%", ], 4),
    Median = round(quantiles["50%", ], 4),
    Q97.5 = round(quantiles["97.5%", ], 4),
    ESS = round(ess, 0),
    ESS_bulk = round(per_parameter[2L, ], 0),
    ESS_tail = round(per_parameter[3L, ], 0),
    MCSE = round(mcse, 4),
    IAT = round(per_parameter[4L, ], 2),
    Rhat = round(per_parameter[5L, ], 3),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

#' Calculate MCMC diagnostics
#'
#' Computes posterior summaries, effective sample sizes, Monte Carlo standard
#' errors (MCSE), integrated autocorrelation times (IAT), and R-hat values
#' after discarding a burn-in fraction.
#'
#' The convergence diagnostics follow Vehtari et al. (2021) and agree with
#' those of the \pkg{posterior} package:
#' * `ESS` is the effective sample size for the mean, estimated from split
#'   chains with Geyer's (1992) initial monotone sequence; `MCSE = SD /
#'   sqrt(ESS)` and `IAT = N / ESS`, where `N` is the number of retained
#'   draws. The IAT is the integrated autocorrelation time used by
#'   Christen & Fox (2010) to compare samplers.
#' * `ESS_bulk` and `ESS_tail` are the rank-normalized bulk and tail
#'   effective sample sizes.
#' * `Rhat` is the rank-normalized split-R-hat (the maximum of its bulk and
#'   tail versions). Because chains are split in halves, it is also defined
#'   for a single chain. Values above 1.01 indicate lack of convergence.
#'
#' When chains have different lengths, ESS and R-hat use the first draws of
#' each chain up to the length of the shortest one.
#'
#' @param samples A matrix of MCMC samples (iterations x parameters), or a
#'   list of matrices for multiple independent chains. For a list, burn-in is
#'   discarded separately from each chain.
#' @param burnin_frac Fraction of samples to discard as burn-in.
#' @param param_names Optional character vector of parameter names. By
#'   default, the column names of `samples` are used, if any.
#' @param title Optional title for printed output.
#'
#' @return A data frame with posterior summaries, ESS (bulk and tail), MCSE,
#'   IAT, and R-hat, returned invisibly after being printed.
#'
#' @references
#' Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Buerkner, P.-C.
#' (2021). Rank-normalization, folding, and localization: An improved R-hat
#' for assessing convergence of MCMC. *Bayesian Analysis*, 16(2), 667--718.
#' \doi{10.1214/20-BA1221}
#'
#' Geyer, C. J. (1992). Practical Markov chain Monte Carlo. *Statistical
#' Science*, 7(4), 473--483.
#'
#' @examples
#' log_post <- function(x) dnorm(x, log = TRUE)
#' res <- twalk(log_post, n_iter = 2000, x0 = -2, xp0 = 2)
#' calculate_diagnostics(
#'   res$samples,
#'   burnin_frac = 0.2,
#'   param_names = "theta",
#'   title = "Standard normal"
#' )
#'
#' @export

calculate_diagnostics <- function(samples, burnin_frac = 0.2, param_names = NULL, title = "") {
  post_burnin_chains <- .prepare_mcmc_chains(samples, burnin_frac)
  results_table <- .diagnostics_table(post_burnin_chains, param_names)

  cat("\nCONVERGENCE DIAGNOSTICS -", title, "\n")
  print(results_table, row.names = FALSE)
  return(invisible(results_table))
}

#' Visualize MCMC results
#'
#' Produces trace plots, marginal densities and joint plots
#' depending on the dimension of the parameter space.
#'
#' @param samples Matrix of MCMC samples.
#' @param true_values Optional vector of true parameter values.
#' @param title Plot title.
#' @param burnin_frac Burn-in fraction.
#' @param true_covariance Optional true covariance matrix.
#' @param show_acf Logical; whether to display autocorrelation plots.
#'
#' @return No return value; called for its side effect of drawing plots on
#'   the current graphics device.
#'
#' @examples
#' log_post <- function(x) dnorm(x, log = TRUE)
#' res <- twalk(log_post, n_iter = 2000, x0 = -2, xp0 = 2)
#' visualize_results(
#'   res$samples,
#'   true_values = 0,
#'   title = "Standard normal"
#' )
#'
#' @export

visualize_results <- function(samples, true_values = NULL,
                              title = "Results", burnin_frac = 0.2,
                              true_covariance = NULL, show_acf = TRUE) {
  post_burnin_chains <- .prepare_mcmc_chains(samples, burnin_frac)
  if (length(post_burnin_chains) != 1L) {
    stop("`visualize_results()` currently accepts one chain at a time.", call. = FALSE)
  }

  analysis_samples <- post_burnin_chains[[1L]]
  n_dim <- ncol(analysis_samples)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))

  colors <- c("#2E86AB", "#A23B72", "#F18F01", "#C73E1D", "#7209B7", "#2F9599")
  bg_color <- "#FAFAFA"

  if (n_dim == 1) {
    layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE))
    par(mar = c(4, 4, 3, 1), bg = bg_color, col.axis = "gray30", col.lab = "gray30", col.main = "gray20")

    hist(analysis_samples[,1], main = paste(title, "- Distribution"),
         xlab = bquote(theta[1]), freq = FALSE, col = colors[1], border = "white",
         breaks = 30, las = 1)

    plot(analysis_samples[, 1], type = 'l', col = colors[1], lwd = 0.8,
         main = paste(title, "- Trace"), xlab = "Iteration", ylab = bquote(theta[1]), las = 1)

    if (show_acf) {
      acf(analysis_samples[,1], main = paste(title, "- Autocorrelation"),
          col = colors[2], lwd = 2, las = 1)
    }

  } else if (n_dim == 2) {
    layout(matrix(c(1, 2, 3, 4), 2, 2, byrow = TRUE))
    par(mar = c(4, 4, 3, 1), bg = bg_color, col.axis = "gray30", col.lab = "gray30", col.main = "gray20")

    plot(analysis_samples, pch = ".", col = adjustcolor(colors[1], alpha.f = 0.6),
         main = title, xlab = bquote(theta[1]), ylab = bquote(theta[2]), las = 1,
         cex = 1.2)
    if (!is.null(true_covariance) && requireNamespace("ellipse", quietly = TRUE)) {
      lines(ellipse::ellipse(true_covariance, centre = true_values), col = colors[2], lwd = 2)
    }

    plot(analysis_samples[, 1], type = 'l', col = colors[1], lwd = 0.8,
         main = bquote(paste("Trace ", theta[1])), xlab = "Iteration", ylab = bquote(theta[1]), las = 1)

    plot(analysis_samples[, 2], type = 'l', col = colors[2], lwd = 0.8,
         main = bquote(paste("Trace ", theta[2])), xlab = "Iteration", ylab = bquote(theta[2]), las = 1)

    dens1 <- density(analysis_samples[,1])
    plot(dens1, main = bquote(paste("Marginal Density ", theta[1])), xlab = bquote(theta[1]),
         col = colors[1], lwd = 3, las = 1, zero.line = FALSE)
    polygon(dens1, col = adjustcolor(colors[1], alpha.f = 0.3), border = NA)

  } else {
    n_plots <- min(6, n_dim)
    if (n_plots <= 4) {
      layout(matrix(1:4, 2, 2, byrow = TRUE))
    } else {
      layout(matrix(1:6, 3, 2, byrow = TRUE))
    }

    par(mar = c(4, 4, 3, 1), bg = bg_color, col.axis = "gray30", col.lab = "gray30", col.main = "gray20")

    for (i in 1:n_plots) {
      if (i <= 3) {
        plot(analysis_samples[, i], type = 'l', col = colors[i], lwd = 0.8,
             main = bquote(paste("Trace ", theta[.(i)])),
             xlab = "Iteration", ylab = bquote(theta[.(i)]), las = 1)
      } else {
        dens <- density(analysis_samples[,i])
        plot(dens, main = bquote(paste("Density ", theta[.(i)])),
             xlab = bquote(theta[.(i)]), col = colors[i], lwd = 3, las = 1, zero.line = FALSE)
        polygon(dens, col = adjustcolor(colors[i], alpha.f = 0.3), border = NA)
      }
    }

    if (n_dim > 6) {
      cat("\nNote: Only showing first 6 dimensions. Total dimensions:", n_dim, "\n")
    }
  }
}
