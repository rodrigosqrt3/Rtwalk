# Less common branches of the drivers and print methods.

test_that("print handles objects without iteration or chain counts", {
  legacy <- structure(list(samples = matrix(0, 5, 2)), class = "twalk")
  expect_output(print(legacy), "Iterations per chain: 5")
})

test_that("multi-chain summaries and errors are reported", {
  log_post <- function(x) -0.5 * sum(x^2)
  set.seed(1)
  multi <- twalk(log_post, n_iter = 200, x0 = rbind(c(-1, 1), c(1, -1)),
                 xp0 = rbind(c(1, -1), c(-1, 1)), n_chains = 2, n_cores = 1,
                 show_progress = FALSE)
  expect_output(print(summary(multi)), "2 chains")

  throwing <- function(x) {
    if (x[1] < -0.5) stop("too negative")
    -0.5 * sum(x^2)
  }
  set.seed(2)
  expect_warning(
    with_errors <- twalk(throwing, n_iter = 300, x0 = rbind(c(0, 0), c(1, 1)),
                         xp0 = rbind(c(1, 1), c(0, 0)), n_chains = 2,
                         n_cores = 1, show_progress = FALSE),
    "too negative"
  )
  expect_identical(with_errors$first_error, "too negative")
})

test_that("multi-chain starting points outside the support are reported", {
  positive <- function(x) if (any(x <= 0)) -Inf else 0
  expect_error(
    twalk(positive, n_iter = 10, x0 = rbind(c(1, 1), c(-1, 1)),
          xp0 = rbind(c(2, 2), c(2, 2)), n_chains = 2, n_cores = 1,
          show_progress = FALSE),
    "Initial points of chain 2"
  )

  failing <- function(x) if (x[1] < 0) stop("bad start") else 0
  expect_error(
    twalk(failing, n_iter = 10, x0 = rbind(c(1, 1), c(-1, 1)),
          xp0 = rbind(c(2, 2), c(2, 2)), n_chains = 2, n_cores = 1,
          show_progress = FALSE),
    "First error: bad start"
  )
})

test_that("unnamed samples keep the log posterior in posterior draws", {
  skip_if_not_installed("posterior")
  unnamed <- structure(
    list(samples = matrix(rnorm(40), ncol = 2), log_posterior = rnorm(20)),
    class = "twalk"
  )
  draws <- posterior::as_draws_array(unnamed)
  expect_identical(posterior::variables(draws), c("theta1", "theta2", "lp__"))
})

test_that("twalk_continue chooses cores automatically for several chains", {
  log_post <- function(x) -0.5 * sum(x^2)
  set.seed(3)
  multi <- twalk(log_post, n_iter = 20, x0 = rbind(c(-1, 1), c(1, -1)),
                 xp0 = rbind(c(1, -1), c(-1, 1)), n_chains = 2, n_cores = 1,
                 show_progress = FALSE)
  more <- twalk_continue(multi, log_post, n_iter = 10, show_progress = FALSE)
  expect_equal(more$n_iter, 30)
})
