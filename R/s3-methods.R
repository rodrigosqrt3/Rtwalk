#' Methods for objects of class 'twalk'
#'
#' S3 methods for inspecting t-walk MCMC output objects.
#'
#' @name twalk-methods
#' @docType methods
#' @keywords internal
NULL

#' @rdname twalk-methods
#' @method print twalk
#' @export
print.twalk <- function(x, ...) {
  cat("t-walk MCMC output\n")
  iterations <- if (!is.null(x$n_iter)) x$n_iter else nrow(x$samples)
  cat("Iterations per chain:", iterations, "\n")
  cat("Dimension:", ncol(x$samples), "\n")
  invisible(x)
}

#' @rdname twalk-methods
#' @method summary twalk
#' @export
summary.twalk <- function(object, burnin_frac = 0.2, ...) {
  samples <- if (!is.null(object$individual_chains) &&
                 length(object$individual_chains) > 1L) {
    lapply(object$individual_chains, function(chain) chain$samples)
  } else {
    object$samples
  }

  calculate_diagnostics(
    samples,
    burnin_frac = burnin_frac
  )
}
