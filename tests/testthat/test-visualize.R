test_that("visualize_results generates plots for 1D, 2D, and high-dimensional samples", {
  pdf(file = tempfile())
  on.exit(dev.off())

  set.seed(123)

  samples_1d <- matrix(rnorm(100), ncol = 1)
  expect_no_error(visualize_results(samples_1d, show_acf = TRUE, title = "1D Test"))
  expect_no_error(visualize_results(samples_1d, show_acf = FALSE, title = "1D Test No ACF"))

  samples_2d <- matrix(rnorm(200), ncol = 2)
  cov_mat <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  expect_no_error(
    visualize_results(samples_2d, true_values = c(0, 0), true_covariance = cov_mat, title = "2D Test")
  )

  samples_3d <- matrix(rnorm(300), ncol = 3)
  expect_no_error(visualize_results(samples_3d, title = "3D Test"))

  samples_5d <- matrix(rnorm(500), ncol = 5)
  expect_no_error(visualize_results(samples_5d, title = "5D Test"))

  samples_7d <- matrix(rnorm(700), ncol = 7)
  expect_output(
    visualize_results(samples_7d, title = "7D Test"),
    "Note: Only showing first 6 dimensions"
  )
})

test_that("visualize_results uses shared sample and burn-in validation", {
  pdf(file = tempfile())
  on.exit(dev.off())

  samples <- matrix(rnorm(200), ncol = 2)

  for (burnin in c(-0.1, 1, 1.2, NA_real_)) {
    expect_error(
      visualize_results(samples, burnin_frac = burnin),
      "`burnin_frac` must be one finite number in \\[0, 1\\)"
    )
  }

  expect_error(
    visualize_results(list(samples, samples)),
    "one chain at a time"
  )
})
