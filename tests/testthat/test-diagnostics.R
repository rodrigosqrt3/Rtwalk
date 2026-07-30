test_that("calculate_diagnostics works with default and custom parameter names", {
  set.seed(123)
  samples <- matrix(rnorm(200), ncol = 2)

  out1 <- expect_output(
    res1 <- calculate_diagnostics(samples, burnin_frac = 0.2, title = "Test Run")
  )
  expect_equal(res1$Parameter, c("theta1", "theta2"))
  expect_equal(nrow(res1), 2)

  out2 <- expect_output(
    res2 <- calculate_diagnostics(samples, burnin_frac = 0.1, param_names = c("alpha", "beta"))
  )
  expect_equal(res2$Parameter, c("alpha", "beta"))
})
