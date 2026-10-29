# Extracted from test-methods.R:18

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "Rtwalk", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
log_post <- function(x) -0.5 * sum(x^2)
res <- twalk(log_posterior = log_post, n_iter = 50, x0 = c(0, 0), xp0 = c(1, 1))
print_output <- capture.output(print(res))
expect_match(print_output, "t-walk MCMC output", all = FALSE)
expect_match(print_output, "Iterations per chain: 50", all = FALSE)
expect_match(print_output, "Chains: 1", all = FALSE)
expect_match(print_output, "Acceptance rate:", all = FALSE)
expect_match(print_output, "Actual move rate:", all = FALSE)
expect_match(print_output, "Kernel statistics: available", all = FALSE)
summary_output <- capture.output(
    sum_res <- summary(res, burnin_frac = 0.2)
  )
expect_match(summary_output, "SAMPLER PERFORMANCE", all = FALSE)
expect_match(summary_output, "Performance by kernel", all = FALSE)
expect_match(summary_output, "traverse", all = FALSE)
