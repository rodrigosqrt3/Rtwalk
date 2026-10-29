# Extracted from test-diagnostics.R:16

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "Rtwalk", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
set.seed(123)
samples <- matrix(rnorm(200), ncol = 2)
out1 <- expect_output(
    res1 <- calculate_diagnostics(samples, burnin_frac = 0.2, title = "Test Run")
  )
expect_equal(res1$Parameter, c("theta1", "theta2"))
expect_equal(nrow(res1), 2)
expect_named(
    res1,
    c(
      "Parameter", "Mean", "SD", "Q2.5", "Median", "Q97.5",
      "ESS", "MCSE", "IAT", "Rhat"
    )
  )
