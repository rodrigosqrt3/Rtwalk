log_std_normal <- function(x) -0.5 * sum(x^2)

test_that("twalk_control validates and normalizes its arguments", {
  control <- twalk_control()
  expect_s3_class(control, "twalk_control")
  expect_equal(unname(control$kernel_probs), c(0.4918, 0.4918, 0.0082, 0.0082))
  expect_equal(c(control$n1, control$aw, control$at), c(4, 1.5, 6))
  expect_output(print(control), "n1 = 4, aw = 1.5, at = 6")

  reordered <- twalk_control(kernel_probs = c(hop = 1, blow = 1, walk = 2, traverse = 4))
  expect_equal(reordered$kernel_probs,
               c(traverse = 0.5, walk = 0.25, blow = 0.125, hop = 0.125))

  expect_error(twalk_control(n1 = 0), "`n1`")
  expect_error(twalk_control(aw = -1), "`aw`")
  expect_error(twalk_control(at = 1), "`at`")
  expect_error(twalk_control(kernel_probs = c(1, 1)), "`kernel_probs`")
  expect_error(twalk_control(kernel_probs = c(a = 1, b = 1, c = 1, d = 1)), "Names")
  expect_error(twalk_control(store_companion = NA), "`store_companion`")
  expect_error(twalk_control(engine = "fortran"))
})

test_that("kernel_probs controls which moves are proposed", {
  set.seed(1)
  walk_only <- twalk(log_std_normal, n_iter = 200, x0 = c(0, 0), xp0 = c(1, 1),
                     show_progress = FALSE,
                     control = twalk_control(kernel_probs = c(0, 1, 0, 0)))
  expect_equal(walk_only$kernel_statistics$proposed, c(0L, 200L, 0L, 0L))
})

test_that("thinning stores every thin-th iteration of the same chain", {
  set.seed(5)
  full <- twalk(log_std_normal, n_iter = 100, x0 = c(0, 0), xp0 = c(1, 1),
                show_progress = FALSE)
  set.seed(5)
  thinned <- twalk(log_std_normal, n_iter = 100, x0 = c(0, 0), xp0 = c(1, 1),
                   thin = 7, show_progress = FALSE)

  keep <- seq(7, 100, by = 7)
  expect_equal(thinned$samples, full$samples[keep, ])
  expect_equal(thinned$companion_samples, full$companion_samples[keep, ])
  expect_equal(thinned$log_posterior, full$log_posterior[keep])
  expect_identical(thinned$kernel_statistics, full$kernel_statistics)
  expect_equal(thinned$n_iter, 100)
  expect_equal(thinned$thin, 7)
})

test_that("the companion trajectory can be dropped to save memory", {
  set.seed(5)
  result <- twalk(log_std_normal, n_iter = 50, x0 = c(0, 0), xp0 = c(1, 1),
                  show_progress = FALSE,
                  control = twalk_control(store_companion = FALSE))
  expect_null(result$companion_samples)
  expect_null(result$companion_log_posterior)
  expect_equal(nrow(result$samples), 50)
  expect_error(coda::as.mcmc(result, companion = TRUE))
})

test_that("twalk_continue extends the same Markov chain", {
  set.seed(11)
  full <- twalk(log_std_normal, n_iter = 300, x0 = c(-1, 1), xp0 = c(1, -1),
                show_progress = FALSE)

  set.seed(11)
  first <- twalk(log_std_normal, n_iter = 100, x0 = c(-1, 1), xp0 = c(1, -1),
                 show_progress = FALSE)
  continued <- twalk_continue(first, log_std_normal, n_iter = 200,
                              show_progress = FALSE)

  expect_equal(continued$samples, full$samples)
  expect_equal(continued$companion_samples, full$companion_samples)
  expect_equal(continued$log_posterior, full$log_posterior)
  expect_equal(continued$kernel_statistics, full$kernel_statistics)
  expect_equal(continued$acceptance_rate, full$acceptance_rate)
  expect_equal(continued$n_iter, 300)

  only_new <- twalk_continue(first, log_std_normal, n_iter = 50,
                             show_progress = FALSE, append = FALSE)
  expect_equal(nrow(only_new$samples), 50)
})

test_that("twalk_continue works for multiple chains and validates input", {
  set.seed(12)
  multi <- twalk(log_std_normal, n_iter = 60,
                 x0 = rbind(c(-2, 2), c(2, -2)), xp0 = rbind(c(-1, 1), c(1, -1)),
                 n_chains = 2, n_cores = 1, show_progress = FALSE)
  more <- twalk_continue(multi, log_std_normal, n_iter = 40, n_cores = 1,
                         show_progress = FALSE)

  expect_equal(more$n_chains, 2)
  expect_equal(nrow(more$individual_chains[[1]]$samples), 100)
  expect_equal(more$individual_chains[[2]]$samples[1:60, ],
               multi$individual_chains[[2]]$samples)
  expect_equal(more$total_iterations, 200)

  expect_error(twalk_continue(list(), log_std_normal, 10), "`twalk` result")
  expect_error(twalk_continue(multi, log_std_normal, 0), "`n_iter`")
  legacy <- structure(list(samples = matrix(0, 2, 1)), class = "twalk")
  expect_error(twalk_continue(legacy, log_std_normal, 10), "older version")
})

test_that("multi-chain results do not depend on the number of cores", {
  run <- function(n_cores) {
    set.seed(99)
    twalk(log_std_normal, n_iter = 40,
          x0 = rbind(c(-2, 2), c(2, -2)), xp0 = rbind(c(-1, 1), c(1, -1)),
          n_chains = 2, n_cores = n_cores, show_progress = FALSE)
  }

  sequential <- run(1)
  parallel_run <- run(2)

  expect_identical(sequential$samples, parallel_run$samples)
  expect_identical(sequential$kernel_statistics, parallel_run$kernel_statistics)
})

test_that("chains use their own starting points and the user's RNG is restored", {
  set.seed(7)
  before_kind <- RNGkind()
  result <- twalk(log_std_normal, n_iter = 1,
                  x0 = rbind(c(-50, 50), c(50, -50)),
                  xp0 = rbind(c(-49, 49), c(49, -49)),
                  n_chains = 2, n_cores = 1, show_progress = FALSE)
  expect_identical(RNGkind(), before_kind)

  first_rows <- t(vapply(result$individual_chains,
                         function(chain) chain$final_state$x, numeric(2)))
  expect_true(all(abs(first_rows) > 40))

  expect_error(
    twalk(log_std_normal, n_iter = 5, x0 = rbind(c(0, 0), c(0, 0)),
          xp0 = rbind(c(1, 1), c(0, 0)), n_chains = 2, show_progress = FALSE),
    "must be different initial points"
  )
  expect_message(
    twalk(log_std_normal, n_iter = 5, x0 = c(0, 0), xp0 = c(1, 1),
          n_chains = 2, n_cores = 1),
    "all chains start from the same points"
  )
})

test_that("starting points outside the support are reported per chain", {
  positive <- function(x) if (any(x <= 0)) -Inf else 0
  expect_error(
    twalk(positive, n_iter = 5, x0 = rbind(c(1, 1), c(-1, 1)),
          xp0 = rbind(c(2, 2), c(2, 2)), n_chains = 2, show_progress = FALSE),
    "chain 2 are outside the support"
  )
})

test_that("an interrupted run returns the completed iterations", {
  calls <- 0L
  interrupting <- function(x) {
    calls <<- calls + 1L
    if (calls == 300L) {
      signalCondition(structure(
        class = c("interrupt", "condition"),
        list(message = "", call = NULL)
      ))
    }
    -0.5 * sum(x^2)
  }

  set.seed(4)
  expect_warning(
    result <- twalk(interrupting, n_iter = 1000, x0 = c(0, 0), xp0 = c(1, 1),
                    show_progress = FALSE),
    "interrupted"
  )
  expect_true(result$interrupted)
  expect_lt(result$n_iter, 1000)
  expect_equal(nrow(result$samples), result$n_iter)
  expect_equal(sum(result$kernel_statistics$proposed), result$n_iter + 1L)
})

test_that("twalk_step performs one transition of the kernel", {
  set.seed(3)
  step <- twalk_step(log_std_normal, x = c(0, 0), xp = c(1, 1))
  expect_named(step, c("x", "xp", "log_x", "log_xp", "kernel", "accepted", "moved", "error"))
  expect_equal(step$log_x, log_std_normal(step$x))
  expect_equal(step$log_xp, log_std_normal(step$xp))
  expect_true(step$kernel %in% c("traverse", "walk", "blow", "hop"))

  # The same seed gives the first transition of twalk()
  set.seed(3)
  run <- twalk(log_std_normal, n_iter = 1, x0 = c(0, 0), xp0 = c(1, 1),
               show_progress = FALSE)
  set.seed(3)
  step <- twalk_step(log_std_normal, x = c(0, 0), xp = c(1, 1),
                     log_x = 0, log_xp = -1)
  expect_equal(step$x, run$samples[1, ])
  expect_equal(step$xp, run$companion_samples[1, ])

  expect_error(twalk_step(log_std_normal, x = 0, xp = c(1, 1)), "same length")
  expect_error(twalk_step(function(x) -Inf, x = 0, xp = 1), "inside the support")
})

test_that("twalk_step can be used inside a Gibbs sampler", {
  set.seed(1)
  y <- rnorm(50, 2, 1.5)
  log_cond_mu <- function(mu, tau) sum(dnorm(y, mu, 1 / sqrt(tau), log = TRUE))

  mu <- 0
  mu_p <- 1
  tau <- 1
  draws <- matrix(NA_real_, 3000, 2)
  for (i in seq_len(3000)) {
    step <- twalk_step(log_cond_mu, mu, mu_p, tau = tau)
    mu <- step$x
    mu_p <- step$xp
    tau <- rgamma(1, 1 + length(y) / 2, 1 + sum((y - mu)^2) / 2)
    draws[i, ] <- c(mu, tau)
  }

  expect_equal(mean(draws[-(1:500), 1]), mean(y), tolerance = 0.1)
})

test_that("posterior draws conversion keeps chains and adds lp__", {
  skip_if_not_installed("posterior")

  set.seed(8)
  single <- twalk(log_std_normal, n_iter = 100, x0 = c(a = 0, b = 0),
                  xp0 = c(1, 1), show_progress = FALSE)
  draws <- posterior::as_draws_array(single, burnin_frac = 0.2)
  expect_s3_class(draws, "draws_array")
  expect_identical(posterior::variables(draws), c("a", "b", "lp__"))
  expect_equal(posterior::niterations(draws), 80)
  expect_equal(
    as.numeric(posterior::extract_variable(draws, "a")),
    unname(single$samples[21:100, "a"])
  )

  multi <- twalk(log_std_normal, n_iter = 50,
                 x0 = rbind(c(0, 0), c(2, 2)), xp0 = rbind(c(1, 1), c(-1, -1)),
                 n_chains = 2, n_cores = 1, show_progress = FALSE)
  draws_df <- posterior::as_draws_df(multi)
  expect_s3_class(draws_df, "draws_df")
  expect_equal(posterior::nchains(draws_df), 2)
  expect_identical(posterior::variables(draws_df), c("theta1", "theta2", "lp__"))
  expect_s3_class(posterior::as_draws(multi, companion = TRUE), "draws_array")
})
