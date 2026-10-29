# Extracted from test-methods.R:67

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "Rtwalk", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
chain1 <- matrix(1:10, ncol = 1)
chain2 <- matrix(101:110, ncol = 1)
object <- structure(
    list(
      samples = rbind(chain1, chain2),
      companion_samples = matrix(numeric(0), ncol = 1),
      all_samples = rbind(chain1, chain2),
      n_iter = 10,
      n_dim = 1,
      individual_chains = list(
        list(samples = chain1),
        list(samples = chain2)
      )
    ),
    class = "twalk"
  )
expect_output(
    result <- summary(object, burnin_frac = 0.2),
    "CONVERGENCE DIAGNOSTICS"
  )
