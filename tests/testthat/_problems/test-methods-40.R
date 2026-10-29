# Extracted from test-methods.R:40

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "Rtwalk", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
primary <- matrix(1:10, ncol = 1)
companion <- matrix(101:110, ncol = 1)
object <- structure(
    list(
      samples = primary,
      companion_samples = companion,
      all_samples = rbind(primary, companion),
      n_iter = 10,
      n_dim = 1
    ),
    class = "twalk"
  )
expect_output(
    result <- summary(object, burnin_frac = 0.2),
    "CONVERGENCE DIAGNOSTICS"
  )
