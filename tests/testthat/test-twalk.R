test_that("twalk runs in sequential mode correctly", {
  log_post <- function(x) -0.5 * sum(x^2)

  set.seed(123)
  res <- twalk(log_posterior = log_post, n_iter = 100, x0 = c(-1, 1), xp0 = c(1, -1), n_chains = 1)

  expect_s3_class(res, "twalk")
  expect_equal(res$n_dim, 2)
  expect_equal(res$n_iter, 100)
  expect_equal(res$n_chains, 1)
  expect_equal(res$total_iterations, 100)
  expect_equal(nrow(res$samples), 100)
  expect_equal(nrow(res$companion_samples), 100)
  expect_equal(nrow(res$all_samples), 200)
  expect_equal(res$all_samples[1:100, ], res$samples)
  expect_equal(res$all_samples[101:200, ], res$companion_samples)
  expect_true(res$acceptance_rate >= 0 && res$acceptance_rate <= 1)
})

test_that("twalk can suppress progress output", {
  log_post <- function(x) -0.5 * sum(x^2)

  expect_silent(
    res <- twalk(
      log_posterior = log_post,
      n_iter = 50,
      x0 = c(0, 0),
      xp0 = c(1, 1),
      n_chains = 1,
      show_progress = FALSE
    )
  )

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

test_that("twalk handles log-posterior errors and rejects invalid outputs", {
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

  expect_error(
    twalk(log_posterior = log_post_vector, n_iter = 10, x0 = c(0, 0), xp0 = c(1, 1)),
    "must return a single numeric value"
  )

  invalid_values <- list(NA_real_, NaN, Inf, "invalid")

  for (value in invalid_values) {
    invalid_log_post <- function(x) value

    expect_error(
      twalk(invalid_log_post, n_iter = 10, x0 = c(0, 0), xp0 = c(1, 1)),
      "must return"
    )
  }
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
  expect_equal(res_par$n_chains, 2)
  expect_equal(res_par$total_iterations, 100)
  expect_equal(nrow(res_par$samples), 100)
  expect_equal(nrow(res_par$companion_samples), 100)
  expect_equal(nrow(res_par$all_samples), 200)
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

test_that("parallel mode respects set.seed", {
  log_post <- function(x) -0.5 * sum(x^2)

  run_parallel <- function() {
    twalk(
      log_posterior = log_post,
      n_iter = 30,
      x0 = c(-1, 1),
      xp0 = c(1, -1),
      n_chains = 2,
      n_cores = 2
    )
  }

  set.seed(123)
  first <- run_parallel()

  set.seed(123)
  second <- run_parallel()

  expect_identical(first$samples, second$samples)
  expect_identical(first$companion_samples, second$companion_samples)
  expect_identical(first$acceptance_rate, second$acceptance_rate)
})

test_that("parallel mode preserves valid initial points inside constrained support", {
  log_post_positive <- function(x) {
    if (any(x <= 0)) {
      return(-Inf)
    }

    sum(stats::dexp(x, rate = 1, log = TRUE))
  }

  set.seed(123)

  expect_no_error(
    result <- twalk(
      log_posterior = log_post_positive,
      n_iter = 20,
      x0 = rep(0.01, 10),
      xp0 = rep(0.02, 10),
      n_chains = 2,
      n_cores = 2
    )
  )

  expect_s3_class(result, "twalk")
})

test_that("parallel mode exports direct log-posterior dependencies", {
  mu_global <- c(0.5, -0.5)
  log_using_global <- function(x) {
    sum(stats::dnorm(x, mean = mu_global, log = TRUE))
  }

  set.seed(123)

  expect_no_error(
    result <- twalk(
      log_posterior = log_using_global,
      n_iter = 20,
      x0 = c(0, 0),
      xp0 = c(1, -1),
      n_chains = 2,
      n_cores = 2,
      show_progress = FALSE
    )
  )

  expect_s3_class(result, "twalk")
})

test_that("parallel mode recursively exports helper dependencies", {
  offset_global <- c(0.5, -0.5)
  helper_global <- function(x) {
    sum(stats::dnorm(x, mean = offset_global, log = TRUE))
  }
  log_using_helper <- function(x) {
    helper_global(x)
  }

  set.seed(123)

  expect_no_error(
    result <- twalk(
      log_posterior = log_using_helper,
      n_iter = 20,
      x0 = c(0, 0),
      xp0 = c(1, -1),
      n_chains = 2,
      n_cores = 2,
      show_progress = FALSE
    )
  )

  expect_s3_class(result, "twalk")
})

test_that("arguments for log_posterior do not alter kernel controls", {
  mu <- 0.9
  n_dim <- 10
  x0 <- rep(0, n_dim)
  xp0 <- rep(1, n_dim)

  log_with_argument <- function(x, p_phi) {
    sum(stats::dnorm(x, mean = p_phi, log = TRUE))
  }

  log_with_closure <- function(x) {
    sum(stats::dnorm(x, mean = mu, log = TRUE))
  }

  set.seed(123)
  with_argument <- twalk(
    log_posterior = log_with_argument,
    n_iter = 50,
    x0 = x0,
    xp0 = xp0,
    p_phi = mu,
    show_progress = FALSE
  )

  set.seed(123)
  with_closure <- twalk(
    log_posterior = log_with_closure,
    n_iter = 50,
    x0 = x0,
    xp0 = xp0,
    show_progress = FALSE
  )

  expect_identical(with_argument$samples, with_closure$samples)
  expect_identical(with_argument$companion_samples, with_closure$companion_samples)
  expect_identical(with_argument$acceptance_rate, with_closure$acceptance_rate)
})

test_that("twalk validates structural inputs", {
  log_post <- function(x) -0.5 * sum(x^2)

  expect_error(
    twalk(log_post, n_iter = 0, x0 = 0, xp0 = 1),
    "`n_iter` must be a positive integer"
  )
  expect_error(
    twalk(log_post, n_iter = -5, x0 = 0, xp0 = 1),
    "`n_iter` must be a positive integer"
  )
  expect_error(
    twalk(log_post, n_iter = 2.5, x0 = 0, xp0 = 1),
    "`n_iter` must be a positive integer"
  )
  expect_error(
    twalk(log_post, n_iter = 10, x0 = c(0, 0), xp0 = c(1, 1, 1)),
    "must have the same length"
  )
  expect_error(
    twalk(log_post, n_iter = 10, x0 = 0, xp0 = 1, n_chains = 0),
    "`n_chains` must be a positive integer"
  )
  expect_error(
    twalk(log_post, n_iter = 10, x0 = 0, xp0 = 1, n_chains = 2, n_cores = 0),
    "`n_cores` must be NULL or a positive integer"
  )
  expect_error(
    twalk(log_post, n_iter = 10, x0 = NA_real_, xp0 = 1),
    "must be non-empty numeric vectors with finite values"
  )
  expect_error(
    twalk(log_post, n_iter = 10, x0 = matrix(0), xp0 = 1),
    "must be non-empty numeric vectors with finite values"
  )
  expect_error(
    twalk(1, n_iter = 10, x0 = 0, xp0 = 1),
    "`log_posterior` must be a function"
  )
  expect_error(
    twalk(log_post, n_iter = 10, x0 = 0, xp0 = 1, show_progress = NA),
    "`show_progress` must be TRUE or FALSE"
  )
})

test_that("internal_call is passed through to log_posterior", {
  log_with_argument <- function(x, internal_call) {
    sum(stats::dnorm(x, mean = internal_call, log = TRUE))
  }
  log_with_closure <- function(x) {
    sum(stats::dnorm(x, mean = 0.5, log = TRUE))
  }

  set.seed(123)
  with_argument <- twalk(
    log_with_argument,
    n_iter = 50,
    x0 = c(0, 0),
    xp0 = c(1, 1),
    show_progress = FALSE,
    internal_call = 0.5
  )

  set.seed(123)
  with_closure <- twalk(
    log_with_closure,
    n_iter = 50,
    x0 = c(0, 0),
    xp0 = c(1, 1),
    show_progress = FALSE
  )

  expect_identical(with_argument$samples, with_closure$samples)
  expect_identical(with_argument$companion_samples, with_closure$companion_samples)
})
