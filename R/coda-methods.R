#' Convert t-walk output to coda objects
#'
#' Converts the primary or companion trajectory from a `twalk` result into
#' the standard MCMC classes supplied by the `coda` package.
#'
#' `as.mcmc()` is intended for single-chain results. For parallel results, use
#' `as.mcmc.list()` so that chain boundaries remain available to convergence
#' diagnostics.
#'
#' @param x An object of class `twalk`.
#' @param burnin_frac Fraction of iterations to discard from the beginning of
#'   each chain.
#' @param companion Logical; if `TRUE`, convert the companion trajectory
#'   instead of the primary trajectory.
#' @param ... Additional arguments reserved for compatibility with the `coda`
#'   generics.
#'
#' @return `as.mcmc()` returns a `coda::mcmc` object. `as.mcmc.list()` returns
#'   a `coda::mcmc.list` object.
#'
#' @examples
#' log_post <- function(x) -0.5 * sum(x^2)
#' set.seed(1)
#' fit <- twalk(log_post, n_iter = 2000, x0 = c(-1, 1), xp0 = c(1, -1),
#'              show_progress = FALSE)
#' chain <- coda::as.mcmc(fit, burnin_frac = 0.2)
#' coda::effectiveSize(chain)
#'
#' @name coda-methods
NULL

.twalk_chains <- function(x, companion = FALSE) {
  if (!is.logical(companion) || length(companion) != 1L || is.na(companion)) {
    stop("`companion` must be TRUE or FALSE.", call. = FALSE)
  }

  trajectory_name <- if (companion) "companion_samples" else "samples"

  if (!is.null(x$individual_chains) && length(x$individual_chains) > 1L) {
    return(lapply(
      x$individual_chains,
      function(chain) chain[[trajectory_name]]
    ))
  }

  list(x[[trajectory_name]])
}

#' @rdname coda-methods
#' @method as.mcmc twalk
#' @importFrom coda as.mcmc
#' @export
as.mcmc.twalk <- function(x, burnin_frac = 0, companion = FALSE, ...) {
  chains <- .twalk_chains(x, companion = companion)

  if (length(chains) != 1L) {
    stop(
      "`as.mcmc()` requires a single-chain result; use `as.mcmc.list()` for parallel chains.",
      call. = FALSE
    )
  }

  prepared_chain <- .prepare_mcmc_chains(chains, burnin_frac)[[1L]]
  coda::mcmc(prepared_chain)
}

#' @rdname coda-methods
#' @method as.mcmc.list twalk
#' @importFrom coda as.mcmc.list
#' @export
as.mcmc.list.twalk <- function(x, burnin_frac = 0, companion = FALSE, ...) {
  chains <- .twalk_chains(x, companion = companion)
  prepared_chains <- .prepare_mcmc_chains(chains, burnin_frac)
  mcmc_chains <- lapply(prepared_chains, coda::mcmc)

  do.call(coda::mcmc.list, mcmc_chains)
}
