# Extracted from test-kernels.R:2

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "Rtwalk", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
expect_equal(dot_product(c(1, 2, 3), c(4, 5, 6)), 32)
