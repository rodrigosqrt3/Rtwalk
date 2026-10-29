# Extracted from test-twalk.R:375

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "Rtwalk", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
log_post <- function(x) -0.5 * sum(x^2)
expect_error(
    twalk(log_post, n_iter = 0, x0 = 0, xp0 = 1),
    "`n_iter` must be a positive integer"
  )
expect_error(
    twalk(log_post, n_iter = -5, x0 = 0, xp0 = 1),
    "`n_iter` must be a positive integer"
  )
expect_error(
    twalk(log_post, n_iter = 2.5, x0 = 0, xp0 = 1),
    "`n_iter` must be a positive integer"
  )
expect_error(
    twalk(log_post, n_iter = 10, x0 = c(0, 0), xp0 = c(1, 1, 1)),
    "must have the same length"
  )
expect_error(
    twalk(log_post, n_iter = 10, x0 = c(0, 0), xp0 = c(0, 0)),
    "must be different initial points"
  )
expect_error(
    twalk(log_post, n_iter = 10, x0 = 0, xp0 = 1, n_chains = 0),
    "`n_chains` must be a positive integer"
  )
expect_error(
    twalk(log_post, n_iter = 10, x0 = 0, xp0 = 1, n_chains = 2, n_cores = 0),
    "`n_cores` must be NULL or a positive integer"
  )
expect_error(
    twalk(log_post, n_iter = 10, x0 = NA_real_, xp0 = 1),
    "must be non-empty numeric vectors with finite values"
  )
expect_error(
    twalk(log_post, n_iter = 10, x0 = matrix(0), xp0 = 1),
    "must be non-empty numeric vectors with finite values"
  )
