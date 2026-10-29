test_that("calculate_diagnostics works with default and custom parameter names", {
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
      "ESS", "ESS_bulk", "ESS_tail", "MCSE", "IAT", "Rhat"
    )
  )
  # Split-R-hat is defined for a single chain
  expect_true(all(is.finite(res1$Rhat)))
  expect_true(all(res1$Rhat < 1.1))

  out2 <- expect_output(
    res2 <- calculate_diagnostics(samples, burnin_frac = 0.1, param_names = c("alpha", "beta"))
  )
  expect_equal(res2$Parameter, c("alpha", "beta"))
})

test_that("calculate_diagnostics reports multi-chain R-hat", {
  set.seed(123)
  converged_chains <- list(
    matrix(rnorm(400), ncol = 2),
    matrix(rnorm(400), ncol = 2)
  )

  expect_output(
    converged <- calculate_diagnostics(converged_chains, burnin_frac = 0),
    "CONVERGENCE DIAGNOSTICS"
  )

  expect_true(all(is.finite(converged$Rhat)))
  expect_true(all(converged$Rhat < 1.1))
  expect_true(all(is.finite(converged$MCSE)))
  expect_true(all(is.finite(converged$IAT)))

  separated_chains <- list(
    matrix(rnorm(200, mean = -3), ncol = 1),
    matrix(rnorm(240, mean = 3), ncol = 1)
  )
  expect_output(
    separated <- calculate_diagnostics(separated_chains, burnin_frac = 0),
    "CONVERGENCE DIAGNOSTICS"
  )

  expect_gt(separated$Rhat, 1.1)
})

test_that("R-hat is unavailable for constant draws", {
  constant_chains <- list(
    matrix(rep(1, 20), ncol = 1),
    matrix(rep(1, 20), ncol = 1)
  )

  expect_true(is.na(.calculate_rhat(constant_chains)))
})

test_that("diagnostics agree with the posterior package", {
  skip_if_not_installed("posterior")

  set.seed(2024)
  ar1 <- function(n, rho, mean = 0) {
    as.numeric(stats::filter(rnorm(n), rho, method = "recursive")) + mean
  }
  scenarios <- list(
    single = matrix(ar1(1001, 0.9), ncol = 1),
    multi = cbind(ar1(500, 0.5), ar1(500, 0.5), ar1(500, 0.7, mean = 0.3)),
    heavy = matrix(rt(800, df = 2), ncol = 2),
    discrete = matrix(rpois(600, 3), ncol = 3)
  )

  for (draws in scenarios) {
    expect_equal(.ess_basic(draws), suppressWarnings(posterior::ess_basic(draws)))
    expect_equal(.ess_bulk(draws), suppressWarnings(posterior::ess_bulk(draws)))
    expect_equal(.ess_tail(draws), suppressWarnings(posterior::ess_tail(draws)))
    expect_equal(.rhat_rank(draws), posterior::rhat(draws))
  }
})

test_that("integrated autocorrelation time matches an AR(1) process", {
  set.seed(1)
  rho <- 0.8
  draws <- matrix(as.numeric(stats::filter(rnorm(2e5), rho, method = "recursive")), ncol = 1)
  iat <- length(draws) / .ess_basic(draws)

  expect_equal(iat, (1 + rho) / (1 - rho), tolerance = 0.05)
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
  expect_error(
    calculate_diagnostics(list(1:10)),
    "numeric matrix or data frame"
  )
  expect_error(
    calculate_diagnostics(matrix(1:2, ncol = 1), burnin_frac = 0.5),
    "leave at least two samples"
  )
})
