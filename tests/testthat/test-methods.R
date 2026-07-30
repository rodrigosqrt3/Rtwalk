test_that("print.twalk and summary.twalk work as expected", {
  log_post <- function(x) -0.5 * sum(x^2)
  res <- twalk(log_posterior = log_post, n_iter = 50, x0 = c(0, 0), xp0 = c(1, 1))

  expect_output(print(res), "t-walk MCMC output")
  expect_output(print(res), "Iterations: 100")

  expect_output(sum_res <- summary(res, burnin_frac = 0.2), "CONVERGENCE DIAGNOSTICS")
  expect_s3_class(sum_res, "data.frame")
})
