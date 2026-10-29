# Extracted from test-engine.R:84

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "Rtwalk", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
run_both_engines <- function(log_post, n_iter, x0, xp0, ...) {
  lapply(c(cpp = "cpp", R = "R"), function(engine) {
    set.seed(20101)
    suppressWarnings(twalk(
      log_post, n_iter = n_iter, x0 = x0, xp0 = xp0, show_progress = FALSE,
      control = twalk_control(engine = engine, ...)
    ))
  })
}
expect_same_chain <- function(runs) {
  expect_equal(runs$cpp$samples, runs$R$samples, tolerance = 1e-10)
  expect_equal(runs$cpp$companion_samples, runs$R$companion_samples, tolerance = 1e-10)
  expect_equal(runs$cpp$log_posterior, runs$R$log_posterior, tolerance = 1e-10)
  expect_identical(runs$cpp$kernel_statistics, runs$R$kernel_statistics)
  expect_identical(runs$cpp$n_evaluations, runs$R$n_evaluations)
  expect_identical(runs$cpp$n_errors, runs$R$n_errors)
}

# test -------------------------------------------------------------------------
seen <- NULL
log_named <- function(theta) {
    seen <<- names(theta)
    dnorm(theta[["mu"]], log = TRUE) + dnorm(theta[["sigma"]], log = TRUE)
  }
set.seed(1)
twalk(log_named, n_iter = 20, x0 = c(mu = 0, sigma = 1), xp0 = c(1, 2),
        show_progress = FALSE)
