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

test_that("calculate_diagnostics validates burn-in and parameter names", {
  samples <- matrix(rnorm(200), ncol = 2)

  for (burnin in c(-0.1, 1, 1.2, NA_real_)) {
    expect_error(
      calculate_diagnostics(samples, burnin_frac = burnin),
      "`burnin_frac` must be one finite number in \\[0, 1\\)"
    )
  }

  expect_error(
    calculate_diagnostics(samples, param_names = "only_one"),
    "one name per parameter"
  )
  expect_error(
    calculate_diagnostics(samples, param_names = c("a", "b", "c")),
    "one name per parameter"
  )
})

test_that("calculate_diagnostics validates chains", {
  expect_error(
    calculate_diagnostics(list()),
    "at least one chain"
  )
  expect_error(
    calculate_diagnostics(list(matrix(1:10, ncol = 1), matrix(1:20, ncol = 2))),
    "same number of columns"
  )
  expect_error(
    calculate_diagnostics(matrix(c(1, NA), ncol = 1)),
    "only finite numeric values"
  )
  expect_error(
    calculate_diagnostics(matrix(1, ncol = 1)),
    "at least two rows"
  )
})
