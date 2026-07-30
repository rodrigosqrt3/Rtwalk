test_that("twalk runs in sequential mode correctly", {
  log_post <- function(x) -0.5 * sum(x^2)

  set.seed(123)
  res <- twalk(log_posterior = log_post, n_iter = 100, x0 = c(-1, 1), xp0 = c(1, -1), n_chains = 1)

  expect_s3_class(res, "twalk")
  expect_equal(res$n_dim, 2)
  expect_equal(res$n_iter, 100)
  expect_equal(nrow(res$all_samples), 200)
  expect_true(res$acceptance_rate >= 0 && res$acceptance_rate <= 1)
})

test_that("twalk handles internal_call flag without progress bar", {
  log_post <- function(x) -0.5 * sum(x^2)

  res <- twalk(log_posterior = log_post, n_iter = 50, x0 = c(0, 0), xp0 = c(1, 1),
               n_chains = 1, internal_call = TRUE)

  expect_s3_class(res, "twalk")
})

test_that("twalk throws error when initial points are outside support", {
  log_post <- function(x) {
    if (x[1] < 0) return(-Inf)
    return(-0.5 * sum(x^2))
  }

  expect_error(
    twalk(log_posterior = log_post, n_iter = 50, x0 = c(-1, 1), xp0 = c(1, 1)),
    "Initial points are outside the support"
  )
})

test_that("twalk handles objective_fun and support_fun edge cases (errors and vector outputs)", {
  log_post_error_at_points <- function(x) {
    if (any(x < 0)) stop("Domain error")
    return(-0.5 * sum(x^2))
  }

  expect_error(
    twalk(log_posterior = log_post_error_at_points, n_iter = 50, x0 = c(-1, -1), xp0 = c(1, 1)),
    "Initial points are outside the support"
  )

  set.seed(42)
  expect_no_error(
    twalk(log_posterior = log_post_error_at_points, n_iter = 50, x0 = c(0.1, 0.1), xp0 = c(1, 1))
  )

  log_post_vector <- function(x) {
    return(c(1, 2))
  }

  expect_no_error(
    twalk(log_posterior = log_post_vector, n_iter = 10, x0 = c(0, 0), xp0 = c(1, 1))
  )
})

test_that("twalk runs in parallel mode correctly", {
  log_post <- function(x) -0.5 * sum(x^2)

  set.seed(123)
  res_par <- twalk(
    log_posterior = log_post,
    n_iter = 50,
    x0 = c(-1, 1),
    xp0 = c(1, -1),
    n_chains = 2,
    n_cores = 2
  )

  expect_s3_class(res_par, "twalk")
  expect_equal(res_par$total_iterations, 100)
  expect_length(res_par$individual_chains, 2)
  expect_true(res_par$acceptance_rate >= 0 && res_par$acceptance_rate <= 1)
})

test_that("twalk parallel mode automatically detects cores when n_cores is NULL", {
  log_post <- function(x) -0.5 * sum(x^2)

  res <- twalk(
    log_posterior = log_post,
    n_iter = 20,
    x0 = c(0, 0),
    xp0 = c(1, 1),
    n_chains = 2,
    n_cores = NULL
  )

  expect_s3_class(res, "twalk")
})
