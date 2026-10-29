#' Methods for objects of class 'twalk'
#'
#' S3 methods for inspecting t-walk MCMC output objects. `print()` displays
#' chain metadata and overall sampler-performance rates. `summary()` computes
#' posterior summaries and convergence diagnostics (see
#' [calculate_diagnostics()]) after discarding a burn-in fraction from each
#' chain, and returns them as a data frame of class `summary.twalk`, whose
#' print method also shows the per-kernel performance of the sampler.
#'
#' @param x,object An object of class `twalk` (or `summary.twalk`).
#' @param burnin_frac Fraction of each chain discarded as burn-in.
#' @param ... Unused.
#'
#' @return `print()` returns `x` invisibly. `summary()` returns a data frame
#'   of class `summary.twalk` with one row per parameter.
#'
#' @examples
#' log_post <- function(x) -0.5 * sum(x^2)
#' set.seed(1)
#' fit <- twalk(log_post, n_iter = 2000, x0 = c(a = -1, b = 1), xp0 = c(1, -1),
#'              show_progress = FALSE)
#' fit
#' s <- summary(fit, burnin_frac = 0.2)
#' s
#' s$Rhat
#'
#' @name twalk-methods
NULL

.print_sampler_performance <- function(object, include_kernels = FALSE) {
  rates <- c(
    acceptance_rate = object$acceptance_rate,
    move_rate = object$move_rate,
    no_move_rate = object$no_move_rate
  )
  has_rates <- length(rates) > 0L
  has_kernel_statistics <- !is.null(object$kernel_statistics)

  if (!has_rates && !(include_kernels && has_kernel_statistics)) {
    return(invisible(NULL))
  }

  cat("\nSAMPLER PERFORMANCE\n")

  if (has_rates) {
    labels <- c(
      acceptance_rate = "Acceptance rate",
      move_rate = "Actual move rate",
      no_move_rate = "Accepted no-move rate"
    )

    for (rate_name in names(rates)) {
      cat(sprintf("%s: %.2f%%\n", labels[[rate_name]], 100 * rates[[rate_name]]))
    }
  }

  if (include_kernels && has_kernel_statistics) {
    kernel_statistics <- object$kernel_statistics
    kernel_statistics$acceptance_rate <- sprintf(
      "%.2f%%",
      100 * kernel_statistics$acceptance_rate
    )
    kernel_statistics$move_rate <- sprintf(
      "%.2f%%",
      100 * kernel_statistics$move_rate
    )

    cat("\nPerformance by kernel:\n")
    print(kernel_statistics, row.names = FALSE)
  }

  invisible(NULL)
}

.twalk_env <- new.env(parent = emptyenv())
.twalk_env$warned_all_samples <- FALSE

#' @export
`$.twalk` <- function(x, name) {
  value <- .subset2(x, name)
  if (!is.null(value) || !identical(name, "all_samples")) {
    return(value)
  }

  # `all_samples` is no longer stored; rebuild it for backward compatibility.
  samples <- .subset2(x, "samples")
  companion <- .subset2(x, "companion_samples")
  if (is.null(samples) || is.null(companion)) {
    return(NULL)
  }
  if (!.twalk_env$warned_all_samples) {
    .twalk_env$warned_all_samples <- TRUE
    warning(
      "`all_samples` is deprecated and will be removed in a future version. ",
      "Use `samples` (the primary trajectory) and `companion_samples` instead. ",
      "This warning is displayed once per session.",
      call. = FALSE
    )
  }
  rbind(samples, companion)
}

#' @rdname twalk-methods
#' @method print twalk
#' @export
print.twalk <- function(x, ...) {
  cat("t-walk MCMC output\n")
  iterations <- if (!is.null(x$n_iter)) x$n_iter else nrow(x$samples)
  chains <- if (!is.null(x$n_chains)) x$n_chains else 1L
  cat("Iterations per chain:", iterations, "\n")
  if (!is.null(x$thin) && x$thin > 1L) {
    cat("Thinning interval:", x$thin, "\n")
  }
  cat("Chains:", chains, "\n")
  cat("Dimension:", ncol(x$samples), "\n")
  .print_sampler_performance(x)

  if (!is.null(x$kernel_statistics)) {
    cat("Kernel statistics: available in `$kernel_statistics`\n")
  }

  invisible(x)
}

.primary_chains <- function(object) {
  if (!is.null(object$individual_chains) &&
      length(object$individual_chains) > 1L) {
    lapply(object$individual_chains, function(chain) chain$samples)
  } else {
    list(object$samples)
  }
}

#' @rdname twalk-methods
#' @method summary twalk
#' @export
summary.twalk <- function(object, burnin_frac = 0.2, ...) {
  chains <- .prepare_mcmc_chains(.primary_chains(object), burnin_frac)
  diagnostics <- .diagnostics_table(chains)

  structure(
    diagnostics,
    class = c("summary.twalk", "data.frame"),
    sampler = list(
      acceptance_rate = object$acceptance_rate,
      move_rate = object$move_rate,
      no_move_rate = object$no_move_rate,
      kernel_statistics = object$kernel_statistics
    ),
    burnin_frac = burnin_frac,
    n_chains = length(chains),
    retained = sum(vapply(chains, nrow, integer(1)))
  )
}

#' @rdname twalk-methods
#' @method print summary.twalk
#' @export
print.summary.twalk <- function(x, ...) {
  sampler <- attr(x, "sampler")
  if (!is.null(sampler)) {
    .print_sampler_performance(sampler, include_kernels = TRUE)
  }

  cat("\nCONVERGENCE DIAGNOSTICS")
  if (!is.null(attr(x, "burnin_frac"))) {
    cat(sprintf(
      " (%d chain%s, %d draws after discarding %g%% burn-in)",
      attr(x, "n_chains"), if (attr(x, "n_chains") > 1L) "s" else "",
      attr(x, "retained"), 100 * attr(x, "burnin_frac")
    ))
  }
  cat("\n")

  table <- x
  class(table) <- "data.frame"
  attributes(table) <- attributes(table)[c("names", "row.names", "class")]
  print(table, row.names = FALSE)

  rhat <- x$Rhat
  if (any(is.finite(rhat) & rhat > 1.01)) {
    cat("\nWarning: R-hat > 1.01 for some parameters; the chains may not have converged.\n")
  }

  invisible(x)
}

#' Plot t-walk output
#'
#' Draws trace plots of the log posterior and of selected parameters, with
#' one colour per chain, next to their marginal densities. The trace of the
#' log posterior is the most useful single summary to assess burn-in.
#'
#' @param x An object of class `twalk`.
#' @param pars Indices or names of the parameters to plot. Defaults to the
#'   first (up to) four parameters.
#' @param burnin_frac Fraction of each chain to omit from the plots.
#' @param ... Unused.
#'
#' @return `x`, invisibly.
#'
#' @examples
#' log_post <- function(x) -0.5 * sum(x^2)
#' fit <- twalk(log_post, n_iter = 2000, x0 = c(a = -1, b = 1),
#'              xp0 = c(1, -1), show_progress = FALSE)
#' plot(fit)
#' @importFrom graphics matplot legend
#' @method plot twalk
#' @export
plot.twalk <- function(x, pars = NULL, burnin_frac = 0, ...) {
  chains <- .prepare_mcmc_chains(.primary_chains(x), burnin_frac)
  n_dim <- ncol(chains[[1L]])
  par_names <- colnames(chains[[1L]])
  if (is.null(par_names)) par_names <- paste0("theta", seq_len(n_dim))

  if (is.null(pars)) {
    pars <- seq_len(min(4L, n_dim))
  } else if (is.character(pars)) {
    pars <- match(pars, par_names)
    if (anyNA(pars)) stop("Unknown parameter names in `pars`.", call. = FALSE)
  } else if (!is.numeric(pars) || any(pars < 1 | pars > n_dim)) {
    stop("`pars` must contain valid parameter indices or names.", call. = FALSE)
  }

  log_post_chains <- if (!is.null(x$individual_chains) &&
                         length(x$individual_chains) > 1L) {
    lapply(x$individual_chains, function(chain) chain$log_posterior)
  } else {
    list(x$log_posterior)
  }
  has_log_post <- all(!vapply(log_post_chains, is.null, logical(1)))

  colors <- c("#2E86AB", "#C73E1D", "#F18F01", "#2F9599", "#7209B7", "#A23B72")
  chain_colors <- rep(colors, length.out = length(chains))

  panels <- c(if (has_log_post) "lp", pars)
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  layout(matrix(seq_len(2L * length(panels)), ncol = 2L, byrow = TRUE),
         widths = c(3, 1))
  par(mar = c(3, 4, 2, 1), mgp = c(2, 0.7, 0), las = 1)

  draw_panel <- function(series, label) {
    n_max <- max(lengths(series))
    padded <- vapply(series, function(s) c(s, rep(NA_real_, n_max - length(s))),
                     numeric(n_max))
    matplot(as.matrix(padded), type = "l", lty = 1, lwd = 0.7,
            col = chain_colors, xlab = "Iteration", ylab = label,
            main = paste("Trace of", label))
    values <- unlist(series)
    values <- values[is.finite(values)]
    dens <- density(values)
    plot(dens, main = "", xlab = label, ylab = "", zero.line = FALSE,
         col = colors[1L], lwd = 2)
    polygon(dens, col = adjustcolor(colors[1L], alpha.f = 0.3), border = NA)
  }

  if (has_log_post) {
    lp <- Map(function(values, chain) {
      n_burn <- length(values) - nrow(chain)
      values[seq.int(n_burn + 1L, length(values))]
    }, log_post_chains, chains)
    draw_panel(lp, "log posterior")
  }
  for (j in pars) {
    draw_panel(lapply(chains, function(chain) chain[, j]), par_names[j])
  }

  invisible(x)
}
