test_that("print.twalk and summary.twalk work as expected", {
  log_post <- function(x) -0.5 * sum(x^2)
  res <- twalk(log_posterior = log_post, n_iter = 50, x0 = c(0, 0), xp0 = c(1, 1))

  expect_output(print(res), "t-walk MCMC output")
  expect_output(print(res), "Iterations per chain: 50")

  expect_output(sum_res <- summary(res, burnin_frac = 0.2), "CONVERGENCE DIAGNOSTICS")
  expect_s3_class(sum_res, "data.frame")
})

test_that("summary uses only the primary time-ordered trajectory", {
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

  expect_equal(result$Mean, 6.5)
})

test_that("summary applies burn-in separately to parallel chains", {
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

  expect_equal(result$Mean, 56.5)
})
