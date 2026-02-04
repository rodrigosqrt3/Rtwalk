#' @importFrom stats sd quantile acf density
#' @importFrom graphics par layout hist lines polygon
#' @importFrom grDevices adjustcolor
NULL

#' Calculate MCMC diagnostics
#'
#' Computes posterior summaries (mean, SD, quantiles) and
#' effective sample sizes after discarding a burn-in fraction.
#'
#' @param samples Matrix of MCMC samples (iterations x parameters).
#' @param burnin_frac Fraction of samples to discard as burn-in.
#' @param param_names Optional character vector of parameter names.
#' @param title Optional title for printed output.
#'
#' @return A data frame with posterior summaries and effective sample sizes.
#'
#' @examples
#' log_post <- function(x) dnorm(x, log = TRUE)
#' res <- twalk(log_post, n_iter = 2000, x0 = -2, xp0 = 2)
#' calculate_diagnostics(
#'   res$all_samples,
#'   burnin_frac = 0.2,
#'   param_names = "theta",
#'   title = "Standard normal"
#' )
#'
#' @export

calculate_diagnostics <- function(samples, burnin_frac = 0.2, param_names = NULL, title = "") {
  n_iter <- nrow(samples)
  n_burnin <- floor(burnin_frac * n_iter)
  post_burnin_samples <- samples[(n_burnin + 1):n_iter, , drop = FALSE]
  n_param <- ncol(samples)

  if (is.null(param_names)) {
    param_names <- paste0("theta", 1:n_param)
  }

  means <- apply(post_burnin_samples, 2, mean)
  sds <- apply(post_burnin_samples, 2, sd)
  quantiles <- apply(post_burnin_samples, 2, quantile, probs = c(0.025, 0.5, 0.975))

  chain <- coda::as.mcmc(post_burnin_samples)
  ess <- coda::effectiveSize(chain)

  results_table <- data.frame(
    Parameter = param_names,
    Mean = round(means, 4),
    SD = round(sds, 4),
    Q2.5 = round(quantiles["2.5%", ], 4),
    Median = round(quantiles["50%", ], 4),
    Q97.5 = round(quantiles["97.5%", ], 4),
    ESS = round(ess, 0)
  )

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
#' @examples
#' log_post <- function(x) dnorm(x, log = TRUE)
#' res <- twalk(log_post, n_iter = 2000, x0 = -2, xp0 = 2)
#' visualize_results(
#'   res$all_samples,
#'   true_values = 0,
#'   title = "Standard normal"
#' )
#'
#' @export

visualize_results <- function(samples, true_values = NULL,
                              title = "Results", burnin_frac = 0.2,
                              true_covariance = NULL, show_acf = TRUE) {
  n_total <- nrow(samples)
  start_index <- ceiling(n_total * burnin_frac) + 1
  analysis_samples <- samples[start_index:n_total, , drop = FALSE]
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

    plot(samples[start_index:n_total, 1], type = 'l', col = colors[1], lwd = 0.8,
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

    plot(samples[start_index:n_total, 1], type = 'l', col = colors[1], lwd = 0.8,
         main = bquote(paste("Trace ", theta[1])), xlab = "Iteration", ylab = bquote(theta[1]), las = 1)

    plot(samples[start_index:n_total, 2], type = 'l', col = colors[2], lwd = 0.8,
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
        plot(samples[start_index:n_total, i], type = 'l', col = colors[i], lwd = 0.8,
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
