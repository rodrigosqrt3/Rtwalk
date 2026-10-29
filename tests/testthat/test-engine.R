# The compiled engine must reproduce the pure-R reference implementation of
# Christen & Fox (2010) draw by draw. On x86-64 the chains are bit-identical;
# a small tolerance allows for fused multiply-add on other platforms.

run_both_engines <- function(log_post, n_iter, x0, xp0, ...) {
  lapply(c(cpp = "cpp", R = "R"), function(engine) {
    set.seed(20101)
    suppressWarnings(twalk(
      log_post, n_iter = n_iter, x0 = x0, xp0 = xp0, show_progress = FALSE,
      control = twalk_control(engine = engine, ...)
    ))
  })
}

expect_same_chain <- function(runs) {
  expect_equal(runs$cpp$samples, runs$R$samples, tolerance = 1e-10)
  expect_equal(runs$cpp$companion_samples, runs$R$companion_samples, tolerance = 1e-10)
  expect_equal(runs$cpp$log_posterior, runs$R$log_posterior, tolerance = 1e-10)
  expect_identical(runs$cpp$kernel_statistics, runs$R$kernel_statistics)
  expect_identical(runs$cpp$n_evaluations, runs$R$n_evaluations)
  expect_identical(runs$cpp$n_errors, runs$R$n_errors)
}

test_that("compiled and reference engines agree on low-dimensional targets", {
  gaussian <- function(x) -0.5 * sum(x^2)
  rosenbrock <- function(x) -(1 / 20) * (100 * (x[2] - x[1]^2)^2 + (1 - x[1])^2)

  expect_same_chain(run_both_engines(gaussian, 500, c(-1, 1), c(1, -1)))
  expect_same_chain(run_both_engines(rosenbrock, 500, c(0, 0), c(1, 1)))
  expect_same_chain(run_both_engines(gaussian, 300, 0.5, 2))
})

test_that("compiled and reference engines agree with partial coordinate updates", {
  scaled <- function(x) -0.5 * sum((x / seq_along(x))^2)
  expect_same_chain(run_both_engines(scaled, 500, rep(0.5, 15), rep(1.5, 15)))
})

test_that("compiled and reference engines agree on constrained supports and errors", {
  positive <- function(x) if (any(x <= 0)) -Inf else sum(dexp(x, log = TRUE))
  throwing <- function(x) {
    if (x[1] < 0) stop("negative")
    -0.5 * sum(x^2)
  }

  expect_same_chain(run_both_engines(positive, 500, rep(0.5, 6), rep(1, 6)))

  runs <- run_both_engines(throwing, 500, c(1, 1), c(2, 2))
  expect_same_chain(runs)
  expect_gt(runs$cpp$n_errors, 0)
})

test_that("compiled and reference engines agree with custom controls", {
  gaussian <- function(x) -0.5 * sum(x^2)
  expect_same_chain(run_both_engines(
    gaussian, 500, rep(0, 6), rep(1, 6),
    n1 = 2, aw = 0.7, at = 3, kernel_probs = c(1, 1, 1, 1)
  ))
})

test_that("the engine evaluates the log posterior once per proposal", {
  evaluations <- 0L
  log_counted <- function(x) {
    evaluations <<- evaluations + 1L
    -0.5 * sum(x^2)
  }

  set.seed(123)
  result <- twalk(log_counted, n_iter = 100, x0 = c(-1, 1), xp0 = c(1, -1),
                  show_progress = FALSE)

  expect_equal(evaluations, result$n_evaluations + 2L)
  expect_lte(result$n_evaluations, result$n_iter)
})

test_that("parameter names are passed to the log posterior", {
  seen <- NULL
  log_named <- function(theta) {
    seen <<- names(theta)
    dnorm(theta[["mu"]], log = TRUE) + dnorm(theta[["sigma"]], log = TRUE)
  }

  set.seed(1)
  twalk(log_named, n_iter = 20, x0 = c(mu = 0, sigma = 1), xp0 = c(1, 2),
        show_progress = FALSE)
  expect_identical(seen, c("mu", "sigma"))
})

test_that("integer and invalid return values are handled", {
  set.seed(1)
  expect_no_error(twalk(function(x) -1L, n_iter = 10, x0 = 0, xp0 = 1,
                        show_progress = FALSE))

  calls <- 0L
  late_failure <- function(x) {
    calls <<- calls + 1L
    if (calls > 5L) NA_real_ else -0.5 * sum(x^2)
  }
  for (engine in c("cpp", "R")) {
    calls <- 0L
    expect_error(
      twalk(late_failure, n_iter = 50, x0 = 0, xp0 = 1, show_progress = FALSE,
            control = twalk_control(engine = engine)),
      "finite value or -Inf"
    )
  }
  expect_error(
    twalk(function(x) factor("a"), n_iter = 5, x0 = 0, xp0 = 1, show_progress = FALSE),
    "single numeric value"
  )
})

test_that("a long run reproduces the moments of a correlated Gaussian", {
  rho <- 0.95
  log_bvn <- function(z) -(z[1]^2 - 2 * rho * z[1] * z[2] + z[2]^2) / (2 * (1 - rho^2))

  set.seed(10)
  result <- twalk(log_bvn, n_iter = 1e5, x0 = c(0, 0), xp0 = c(1, 1),
                  show_progress = FALSE)
  draws <- result$samples[-(1:1e4), ]

  expect_equal(colMeans(draws), c(0, 0), tolerance = 0.15)
  expect_equal(apply(draws, 2, var), c(1, 1), tolerance = 0.1)
  expect_equal(cor(draws)[1, 2], rho, tolerance = 0.01)
})
