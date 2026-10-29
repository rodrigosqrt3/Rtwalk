test_that("print.twalk and summary.twalk work as expected", {
  log_post <- function(x) -0.5 * sum(x^2)
  res <- twalk(log_posterior = log_post, n_iter = 50, x0 = c(0, 0), xp0 = c(1, 1),
               show_progress = FALSE)

  print_output <- capture.output(print(res))
  expect_match(print_output, "t-walk MCMC output", all = FALSE)
  expect_match(print_output, "Iterations per chain: 50", all = FALSE)
  expect_match(print_output, "Chains: 1", all = FALSE)
  expect_match(print_output, "Acceptance rate:", all = FALSE)
  expect_match(print_output, "Actual move rate:", all = FALSE)
  expect_match(print_output, "Kernel statistics: available", all = FALSE)

  sum_res <- summary(res, burnin_frac = 0.2)
  expect_s3_class(sum_res, "summary.twalk")
  expect_s3_class(sum_res, "data.frame")
  expect_equal(nrow(sum_res), 2)

  summary_output <- capture.output(print(sum_res))
  expect_match(summary_output, "SAMPLER PERFORMANCE", all = FALSE)
  expect_match(summary_output, "Performance by kernel", all = FALSE)
  expect_match(summary_output, "traverse", all = FALSE)
  expect_match(summary_output, "CONVERGENCE DIAGNOSTICS", all = FALSE)
  expect_match(summary_output, "40 draws after discarding 20% burn-in", all = FALSE)
})

test_that("print reports thinning and summary flags unconverged chains", {
  log_post <- function(x) -0.5 * sum(x^2)
  set.seed(1)
  res <- twalk(log_post, n_iter = 100, x0 = c(0, 0), xp0 = c(1, 1),
               thin = 5, show_progress = FALSE)
  expect_match(capture.output(print(res)), "Thinning interval: 5", all = FALSE)

  stuck <- structure(
    list(samples = matrix(c(1:50, rnorm(50)), ncol = 1)),
    class = "twalk"
  )
  expect_output(print(summary(stuck, burnin_frac = 0)), "R-hat > 1.01")
})

test_that("summary uses only the primary time-ordered trajectory", {
  primary <- matrix(1:10, ncol = 1)
  companion <- matrix(101:110, ncol = 1)
  object <- structure(
    list(
      samples = primary,
      companion_samples = companion,
      n_iter = 10,
      n_dim = 1
    ),
    class = "twalk"
  )

  result <- summary(object, burnin_frac = 0.2)
  expect_output(print(result), "CONVERGENCE DIAGNOSTICS")
  expect_equal(result$Mean, 6.5)
})

test_that("summary applies burn-in separately to parallel chains", {
  chain1 <- matrix(1:10, ncol = 1)
  chain2 <- matrix(101:110, ncol = 1)

  object <- structure(
    list(
      samples = rbind(chain1, chain2),
      companion_samples = matrix(numeric(0), ncol = 1),
      n_iter = 10,
      n_dim = 1,
      individual_chains = list(
        list(samples = chain1),
        list(samples = chain2)
      )
    ),
    class = "twalk"
  )

  result <- summary(object, burnin_frac = 0.2)
  expect_equal(result$Mean, 56.5)
})

test_that("summary uses parameter names from the starting point", {
  log_post <- function(theta) -0.5 * sum(theta^2)
  set.seed(2)
  res <- twalk(log_post, n_iter = 200, x0 = c(alpha = -1, beta = 1),
               xp0 = c(1, -1), show_progress = FALSE)

  expect_identical(colnames(res$samples), c("alpha", "beta"))
  expect_identical(summary(res)$Parameter, c("alpha", "beta"))
})

test_that("all_samples is computed on demand with a deprecation warning", {
  log_post <- function(x) -0.5 * sum(x^2)
  set.seed(3)
  res <- twalk(log_post, n_iter = 30, x0 = c(0, 0), xp0 = c(1, 1),
               show_progress = FALSE)

  expect_null(.subset2(res, "all_samples"))

  .twalk_env$warned_all_samples <- FALSE
  expect_warning(all_samples <- res$all_samples, "deprecated")
  expect_equal(all_samples, rbind(res$samples, res$companion_samples))
  expect_silent(res$all_samples)
})

test_that("plot.twalk draws traces for single and multiple chains", {
  pdf(file = tempfile())
  on.exit(dev.off())

  log_post <- function(x) -0.5 * sum(x^2)
  set.seed(4)
  single <- twalk(log_post, n_iter = 100, x0 = c(a = 0, b = 0, c = 0),
                  xp0 = c(1, 1, 1), show_progress = FALSE)
  expect_invisible(plot(single))
  expect_no_error(plot(single, pars = "b", burnin_frac = 0.5))
  expect_error(plot(single, pars = "zzz"), "Unknown parameter")
  expect_error(plot(single, pars = 10), "valid parameter")

  multi <- twalk(log_post, n_iter = 100,
                 x0 = rbind(c(0, 0), c(2, 2)), xp0 = rbind(c(1, 1), c(-1, -1)),
                 n_chains = 2, n_cores = 1, show_progress = FALSE)
  expect_no_error(plot(multi))
})
