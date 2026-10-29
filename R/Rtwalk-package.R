#' Rtwalk: An MCMC Sampler Using the t-Walk Algorithm
#'
#' Implements the t-walk of Christen & Fox (2010), a general-purpose MCMC
#' sampler for continuous distributions that requires no tuning. The main
#' function is [twalk()]; see also [twalk_control()], [twalk_continue()],
#' [twalk_step()], [summary.twalk()] and [plot.twalk()].
#'
#' Vignettes:
#' * `vignette("validation_study", package = "Rtwalk")`: correctness checks
#'   against known targets.
#' * `vignette("reproducing_christen_fox", package = "Rtwalk")`: the
#'   experiments of the original article.
#' * `vignette("sampler_comparison", package = "Rtwalk")`: a benchmark
#'   against other general-purpose samplers.
#'
#' @references Christen, J. A. and Fox, C. (2010). A general purpose sampling
#'   algorithm for continuous distributions (the t-walk). *Bayesian Analysis*,
#'   5(2), 263--282. \doi{10.1214/10-BA603}
#'
#' @keywords internal
"_PACKAGE"
