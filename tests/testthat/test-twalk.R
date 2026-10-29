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
  expect_length(res$log_posterior, 100)
  expect_length(res$companion_log_posterior, 100)
  expect_equal(res$log_posterior, apply(res$samples, 1, log_post))
  expect_equal(res$companion_log_posterior, apply(res$companion_samples, 1, log_post))
  expect_true(res$acceptance_rate >= 0 && res$acceptance_rate <= 1)
  expect_true(res$move_rate >= 0 && res$move_rate <= res$acceptance_rate)
  expect_equal(res$no_move_rate, res$acceptance_rate - res$move_rate)
  expect_named(
    res$kernel_statistics,
    c(
      "kernel", "proposed", "accepted", "moved",
      "acceptance_rate", "move_rate"
    )
  )
  expect_equal(sum(res$kernel_statistics$proposed), res$n_iter)
  expect_equal(
    sum(res$kernel_statistics$accepted),
    res$acceptance_rate * res$n_iter
  )
  expect_equal(
    sum(res$kernel_statistics$moved),
    res$move_rate * res$n_iter
  )
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
  expect_warning(
    result <- twalk(
      log_posterior = log_post_error_at_points, n_iter = 200,
      x0 = c(0.1, 0.1), xp0 = c(1, 1), show_progress = FALSE
    ),
    "raised an error in .* First error: Domain error"
  )
  expect_gt(result$n_errors, 0)
  expect_identical(result$first_error, "Domain error")
  expect_true(all(result$samples >= 0))
  expect_equal(sum(result$kernel_statistics$proposed), 200)

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
  expect_length(res_par$log_posterior, 100)
  expect_length(res_par$individual_chains, 2)
  expect_true(res_par$acceptance_rate >= 0 && res_par$acceptance_rate <= 1)
  expect_true(res_par$move_rate >= 0 && res_par$move_rate <= res_par$acceptance_rate)
  expect_equal(res_par$no_move_rate, res_par$acceptance_rate - res_par$move_rate)
  expect_equal(
    sum(res_par$kernel_statistics$proposed),
    res_par$total_iterations
  )
  expect_equal(
    sum(res_par$kernel_statistics$accepted),
    res_par$acceptance_rate * res_par$total_iterations
  )
  expect_equal(
    sum(res_par$kernel_statistics$moved),
    res_par$move_rate * res_par$total_iterations
  )
  expect_true(all(vapply(
    res_par$individual_chains,
    function(chain) {
      !is.null(chain$move_rate) &&
        !is.null(chain$no_move_rate) &&
        !is.null(chain$kernel_statistics) &&
        isTRUE(all.equal(
          chain$move_rate,
          sum(chain$kernel_statistics$moved) / chain$n_iter
        )) &&
        isTRUE(all.equal(
          chain$no_move_rate,
          chain$acceptance_rate - chain$move_rate
        ))
    },
    logical(1)
  )))
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

test_that("core detection falls back safely when detection is unavailable", {
  unavailable_detector <- function() NA_integer_

  expect_equal(.detect_available_cores(unavailable_detector), 1L)
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
  expect_identical(first$move_rate, second$move_rate)
  expect_identical(first$no_move_rate, second$no_move_rate)
  expect_identical(first$kernel_statistics, second$kernel_statistics)
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

test_that("parallel mode accepts functions without an enclosing environment", {
  set.seed(123)

  expect_no_error(
    result <- twalk(
      log_posterior = base::sum,
      n_iter = 10,
      x0 = 0,
      xp0 = 1,
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
    twalk(log_post, n_iter = 10, x0 = c(0, 0), xp0 = c(0, 0)),
    "must be different initial points"
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
    twalk(log_post, n_iter = 10, x0 = matrix(0, 2, 1), xp0 = 1),
    "one row per chain"
  )
  expect_error(
    twalk(log_post, n_iter = 10, x0 = "a", xp0 = 1),
    "must be non-empty numeric vectors with finite values"
  )
  expect_error(
    twalk(log_post, n_iter = 10, x0 = 0, xp0 = 1, thin = 20),
    "`thin` must be a positive integer"
  )
  expect_error(
    twalk(log_post, n_iter = 10, x0 = 0, xp0 = 1, control = list()),
    "`control` must be created with `twalk_control\\(\\)`"
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

test_that("log_posterior is evaluated once per proposal", {
  evaluations <- 0L
  log_counted <- function(x) {
    evaluations <<- evaluations + 1L
    -0.5 * sum(x^2)
  }

  set.seed(123)
  result <- twalk(
    log_counted,
    n_iter = 100,
    x0 = c(-1, 1),
    xp0 = c(1, -1),
    show_progress = FALSE
  )

  expect_s3_class(result, "twalk")
  expect_equal(evaluations, result$n_iter + 2L)
})

test_that("twalk is invariant under scalar affine transformations", {
  log_target <- function(x) {
    mu <- c(0.5, -1, 2)
    sd <- c(1, 2, 0.7)
    -0.5 * sum(((x - mu) / sd)^2)
  }

  scale <- 4.25
  shift <- c(10, -3, 7)
  log_target_transformed <- function(z) {
    log_target((z - shift) / scale) -
      length(z) * log(abs(scale))
  }

  x0 <- c(-1, 0, 1)
  xp0 <- c(2, -2, 3)

  set.seed(20260829)
  original <- twalk(
    log_target,
    n_iter = 500,
    x0 = x0,
    xp0 = xp0,
    show_progress = FALSE
  )

  set.seed(20260829)
  transformed <- twalk(
    log_target_transformed,
    n_iter = 500,
    x0 = scale * x0 + shift,
    xp0 = scale * xp0 + shift,
    show_progress = FALSE
  )

  expected_samples <- sweep(scale * original$samples, 2, shift, "+")
  expected_companion <- sweep(
    scale * original$companion_samples,
    2,
    shift,
    "+"
  )

  expect_equal(transformed$samples, expected_samples, tolerance = 1e-10)
  expect_equal(
    transformed$companion_samples,
    expected_companion,
    tolerance = 1e-10
  )
  expect_identical(transformed$acceptance_rate, original$acceptance_rate)
  expect_identical(transformed$move_rate, original$move_rate)
  expect_identical(transformed$no_move_rate, original$no_move_rate)
})
