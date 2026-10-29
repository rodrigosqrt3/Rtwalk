#' Convert t-walk output to posterior draws objects
#'
#' Methods for the generics of the \pkg{posterior} package, which give access
#' to its summaries, diagnostics and plotting ecosystem (for example
#' `posterior::summarise_draws()` or the \pkg{bayesplot} package). Chain
#' boundaries are preserved, and the log posterior of the stored points is
#' included as the variable `lp__`, following the convention of Stan.
#'
#' These methods are registered when \pkg{posterior} is loaded; the package
#' is only suggested.
#'
#' @param x An object of class `twalk`.
#' @param burnin_frac Fraction of iterations to discard from the beginning of
#'   each chain.
#' @param companion Logical; if `TRUE`, convert the companion trajectory
#'   instead of the primary one.
#' @param ... Unused.
#'
#' @return A `draws_array` (or `draws_df`) object.
#'
#' @examples
#' if (requireNamespace("posterior", quietly = TRUE)) {
#'   log_post <- function(x) -0.5 * sum(x^2)
#'   fit <- twalk(log_post, n_iter = 2000, x0 = c(a = -1, b = 1),
#'                xp0 = c(1, -1), show_progress = FALSE)
#'   draws <- posterior::as_draws_array(fit, burnin_frac = 0.2)
#'   posterior::summarise_draws(draws)
#' }
#' @name posterior-methods
NULL

.twalk_draws_array <- function(x, burnin_frac = 0, companion = FALSE) {
  if (!requireNamespace("posterior", quietly = TRUE)) {
    stop("The 'posterior' package is required for this conversion.", call. = FALSE) # nocov
  }

  chains <- .twalk_chains(x, companion = companion)
  log_post_field <- if (companion) "companion_log_posterior" else "log_posterior"
  log_post <- if (!is.null(x$individual_chains) && length(x$individual_chains) > 1L) {
    lapply(x$individual_chains, function(chain) chain[[log_post_field]])
  } else {
    list(x[[log_post_field]])
  }
  has_log_post <- all(!vapply(log_post, is.null, logical(1)))

  if (has_log_post) {
    chains <- Map(cbind, chains, lp__ = log_post)
  }
  chains <- .prepare_mcmc_chains(chains, burnin_frac)

  n_param <- ncol(chains[[1L]])
  variables <- colnames(chains[[1L]])
  default_names <- paste0("theta", seq_len(n_param - has_log_post))
  if (is.null(variables)) {
    variables <- c(default_names, if (has_log_post) "lp__")
  } else {
    missing_names <- is.na(variables) | variables == ""
    variables[missing_names] <- default_names[missing_names]
  }

  common_length <- min(vapply(chains, nrow, integer(1)))
  values <- array(
    NA_real_,
    dim = c(common_length, length(chains), n_param),
    dimnames = list(NULL, NULL, variables)
  )
  for (i in seq_along(chains)) {
    values[, i, ] <- chains[[i]][seq_len(common_length), , drop = FALSE]
  }

  posterior::as_draws_array(values)
}

#' @rdname posterior-methods
#' @exportS3Method posterior::as_draws_array
as_draws_array.twalk <- function(x, burnin_frac = 0, companion = FALSE, ...) {
  .twalk_draws_array(x, burnin_frac = burnin_frac, companion = companion)
}

#' @rdname posterior-methods
#' @exportS3Method posterior::as_draws
as_draws.twalk <- function(x, burnin_frac = 0, companion = FALSE, ...) {
  .twalk_draws_array(x, burnin_frac = burnin_frac, companion = companion)
}

#' @rdname posterior-methods
#' @exportS3Method posterior::as_draws_df
as_draws_df.twalk <- function(x, burnin_frac = 0, companion = FALSE, ...) {
  posterior::as_draws_df(
    .twalk_draws_array(x, burnin_frac = burnin_frac, companion = companion)
  )
}
