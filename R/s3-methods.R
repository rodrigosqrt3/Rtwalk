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
  cat("Iterations:", nrow(x$all_samples), "\n")
  cat("Dimension:", ncol(x$all_samples), "\n")
  invisible(x)
}

#' @rdname twalk-methods
#' @method summary twalk
#' @export
summary.twalk <- function(object, burnin_frac = 0.2, ...) {
  calculate_diagnostics(
    object$all_samples,
    burnin_frac = burnin_frac
  )
}
