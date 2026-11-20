# Penalised t-walk (Medina-Aguayo & Christen, 2022).

log_std_normal <- function(x) -0.5 * sum(x^2)

# Two-component Gaussian mixture of Example 1 of Medina-Aguayo & Christen
# (2022), with weights `w` and 1 - `w`.
log_mixture <- function(w = 0.5) {
  mu2 <- c(20, -20)
  prec1 <- solve(matrix(c(1, 0.1, 0.1, 1), 2))
  prec2 <- solve(matrix(c(16, 16, 16, 25), 2))
  log_det1 <- -0.5 * log(det(matrix(c(1, 0.1, 0.1, 1), 2)))
  log_det2 <- -0.5 * log(det(matrix(c(16, 16, 16, 25), 2)))
  function(x) {
    d2 <- x - mu2
    a <- log(w) + log_det1 - 0.5 * sum(x * (prec1 %*% x))
    b <- log(1 - w) + log_det2 - 0.5 * sum(d2 * (prec2 %*% d2))
    m <- max(a, b)
    m + log(exp(a - m) + exp(b - m))
  }
}

test_that("without the penalised move the chains equal those of version 2.2.0", {
  reference <- readRDS(test_path("fixtures", "twalk_2.2.0_reference.rds"))
  log_post <- function(x) {
    if (x[1] < -3) return(-Inf)
    -0.5 * sum(x^2 / c(1, 4, 0.25)) - 0.1 * x[1] * x[2]
  }

  for (engine in c("cpp", "R")) {
    control <- twalk_control(engine = engine)
    set.seed(123)
    single <- twalk(log_post, n_iter = 2000, x0 = c(0, 1, -1), xp0 = c(1, -1, 0.5),
                    show_progress = FALSE, control = control)
    set.seed(456)
    multi <- twalk(log_post, n_iter = 1000, x0 = c(0, 1, -1), xp0 = c(1, -1, 0.5),
                   n_chains = 2, thin = 2, show_progress = FALSE, control = control)

    # The reference was generated on x86_64 Windows, where the chains must be
    # bit for bit identical; elsewhere (including Windows on ARM, where clang
    # fuses multiply-adds) the last bits may be rounded differently, so a
    # tolerance is used.
    compare <- if (identical(R.version$platform, "x86_64-w64-mingw32")) {
      expect_identical
    } else {
      function(object, expected) expect_equal(object, expected, tolerance = 1e-10)
    }
    for (field in names(reference$single)) {
      compare(single[[field]], reference$single[[field]])
      compare(multi[[field]], reference$multi[[field]])
    }
  }
})

run_both_engines_penalty <- function(log_post, n_iter, x0, xp0, ...) {
  lapply(c(cpp = "cpp", R = "R"), function(engine) {
    set.seed(2022)
    suppressWarnings(twalk(
      log_post, n_iter = n_iter, x0 = x0, xp0 = xp0, show_progress = FALSE,
      control = twalk_control(engine = engine, ...)
    ))
  })
}

# As in test-engine.R, the floating-point values may differ in the last bits
# between the engines on some platforms (the C++ compiler and R's arithmetic
# need not round alike), so they are compared with a tolerance; the accept and
# reject decisions must still agree exactly.
expect_same_penalised_chain <- function(runs) {
  expect_equal(runs$cpp$samples, runs$R$samples, tolerance = 1e-10)
  expect_equal(runs$cpp$companion_samples, runs$R$companion_samples, tolerance = 1e-10)
  expect_equal(runs$cpp$log_posterior, runs$R$log_posterior, tolerance = 1e-10)
  expect_identical(runs$cpp$kernel_statistics, runs$R$kernel_statistics)
  expect_identical(runs$cpp$n_evaluations, runs$R$n_evaluations)
  expect_identical(runs$cpp$n_errors, runs$R$n_errors)
  expect_gt(runs$cpp$kernel_statistics$proposed[5], 0)
}

test_that("compiled and reference engines agree with the penalised move", {
  positive <- function(x) if (any(x <= 0)) -Inf else sum(dexp(x, log = TRUE))
  throwing <- function(x) {
    if (x[1] < -1) stop("too negative")
    -0.5 * sum(x^2)
  }

  expect_same_penalised_chain(run_both_engines_penalty(
    log_std_normal, 1000, c(-1, 1), c(1, -1), penalty_prob = 0.3
  ))
  expect_same_penalised_chain(run_both_engines_penalty(
    log_mixture(), 1000, c(0, 0), c(1, 1), penalty_prob = 0.5, penalty_kappa = 2
  ))
  expect_same_penalised_chain(run_both_engines_penalty(
    log_std_normal, 500, rep(0.5, 15), rep(1.5, 15), penalty_prob = 0.2
  ))
  expect_same_penalised_chain(run_both_engines_penalty(
    positive, 1000, rep(0.5, 4), rep(1, 4), penalty_prob = 0.5
  ))

  runs <- run_both_engines_penalty(throwing, 1000, c(0, 0), c(1, 1),
                                   penalty_prob = 0.5)
  expect_same_penalised_chain(runs)
  expect_gt(runs$cpp$n_errors, 0)
})

test_that("the penalised move reports its own statistics", {
  set.seed(3)
  fit <- twalk(log_std_normal, n_iter = 2000, x0 = c(-1, 1), xp0 = c(1, -1),
               show_progress = FALSE, control = twalk_control(penalty_prob = 0.2))
  stats <- fit$kernel_statistics

  expect_identical(stats$kernel, c("traverse", "walk", "blow", "hop", "penalty"))
  expect_equal(sum(stats$proposed), 2000)
  expect_equal(stats$proposed[5] / 2000, 0.2, tolerance = 0.15)
  expect_lte(stats$moved[5], stats$accepted[5])
  expect_output(print(summary(fit)), "penalty")

  # Each penalised move evaluates the log posterior at most twice.
  expect_lte(fit$n_evaluations, sum(stats$proposed) + stats$proposed[5])
})

test_that("multiple chains, continuation and thinning work with the move", {
  control <- twalk_control(penalty_prob = 0.2)
  set.seed(21)
  full <- twalk(log_std_normal, n_iter = 300, x0 = c(-1, 1), xp0 = c(1, -1),
                show_progress = FALSE, control = control)
  set.seed(21)
  first <- twalk(log_std_normal, n_iter = 100, x0 = c(-1, 1), xp0 = c(1, -1),
                 show_progress = FALSE, control = control)
  continued <- twalk_continue(first, log_std_normal, n_iter = 200,
                              show_progress = FALSE)
  expect_equal(continued$samples, full$samples)
  expect_equal(continued$kernel_statistics, full$kernel_statistics)

  set.seed(22)
  multi <- twalk(log_std_normal, n_iter = 100, x0 = c(-1, 1), xp0 = c(1, -1),
                 n_chains = 2, thin = 5, show_progress = FALSE, control = control)
  expect_equal(nrow(multi$kernel_statistics), 5)
  expect_equal(sum(multi$kernel_statistics$proposed), 200)
  expect_equal(nrow(multi$samples), 40)
})

test_that("the rejection sampler reproduces Table 1 of the paper", {
  # Acceptance rates of Algorithm 2 for the multivariate t penalty with 2
  # degrees of freedom and a t proposal with 1 degree of freedom.
  table1 <- rbind(
    c(d = 2, kappa = 2, rate = 0.8671),
    c(d = 2, kappa = 3, rate = 0.9275),
    c(d = 2, kappa = 4, rate = 0.9551),
    c(d = 4, kappa = 2, rate = 0.9802),
    c(d = 4, kappa = 3, rate = 0.9931)
  )

  set.seed(4)
  for (i in seq_len(nrow(table1))) {
    d <- table1[i, "d"]
    trials <- vapply(seq_len(20000), function(k) {
      penalty_shift(rep(0, d), rep(1, d), table1[i, "kappa"])$n_trials
    }, integer(1))
    expect_equal(20000 / sum(trials), unname(table1[i, "rate"]), tolerance = 0.005)
  }
})

test_that("the penalised shift follows the penalised density", {
  # In one dimension with scale 1 the shift has density
  # g(t) phi(t) / Z, with g the Cauchy density with scale kappa and
  # phi(t) = 1 - (1 + t^2 / 2)^(-3 / 2).
  # The CDF is the Cauchy CDF minus the integral of the penalty term, which
  # is symmetric and decays fast.
  kappa <- 3
  removed <- function(t) stats::dcauchy(t, scale = kappa) * (1 + t^2 / 2)^(-3 / 2)
  density_unnorm <- function(t) stats::dcauchy(t, scale = kappa) - removed(t)
  half <- stats::integrate(removed, 0, Inf)$value
  z_const <- 1 - 2 * half
  cdf <- function(q) {
    vapply(q, function(v) {
      lower_tail <- if (v < 0) {
        stats::integrate(removed, -v, Inf)$value
      } else {
        half + stats::integrate(removed, 0, v)$value
      }
      (stats::pcauchy(v, scale = kappa) - lower_tail) / z_const
    }, numeric(1))
  }

  set.seed(5)
  shifts <- vapply(seq_len(3000), function(k) penalty_shift(0, 1, kappa)$delta,
                   numeric(1))
  expect_gt(stats::ks.test(shifts, cdf)$p.value, 0.001)
  expect_equal(mean(abs(shifts) < 0.5), stats::integrate(density_unnorm, -0.5, 0.5)$value / z_const,
               tolerance = 0.2)
})

test_that("coordinates with equal values are not moved", {
  set.seed(6)
  shift <- penalty_shift(c(1, 2, 3), c(1, 5, 3), kappa = 3)
  expect_identical(shift$delta[c(1, 3)], c(0, 0))
  expect_true(shift$delta[2] != 0)

  same <- penalty_shift(c(1, 2), c(1, 2))
  expect_identical(same, list(delta = c(0, 0), n_trials = 0L))

  # x = x': the penalised move is the identity, accepted without evaluation.
  for (engine in c("cpp", "R")) {
    evaluations <- 0L
    counted <- function(x) {
      evaluations <<- evaluations + 1L
      -0.5 * sum(x^2)
    }
    set.seed(1)
    step <- twalk_step(counted, c(1, 2), c(1, 2), log_x = -2.5, log_xp = -2.5,
                       control = twalk_control(penalty_prob = 0.999, engine = engine))
    expect_identical(step$kernel, "penalty")
    expect_true(step$accepted)
    expect_false(step$moved)
    expect_identical(evaluations, 0L)
  }
})

test_that("non-finite proposals are rejected without evaluation", {
  for (engine in c("cpp", "R")) {
    set.seed(1)
    step <- twalk_step(function(x) stop("should not be called"),
                       1e308, -1e308, log_x = 0, log_xp = 0,
                       control = twalk_control(penalty_prob = 0.999, engine = engine))
    expect_identical(step$kernel, "penalty")
    expect_false(step$accepted)
    expect_identical(step$x, 1e308)
  }
})

test_that("the penalised t-walk preserves the target", {
  # Mostly penalised moves on a correlated Gaussian.
  log_corr <- function(x) -0.5 * (x[1]^2 + (x[2] - x[1])^2 / 0.25)
  set.seed(7)
  fit <- twalk(log_corr, n_iter = 1e5, x0 = c(0, 0), xp0 = c(1, 1),
               show_progress = FALSE, control = twalk_control(penalty_prob = 0.9))
  draws <- fit$samples[-(1:5000), ]
  expect_equal(colMeans(draws), c(0, 0), tolerance = 0.1)
  expect_equal(c(stats::var(draws)), c(1, 1, 1, 1.25), tolerance = 0.1)

  # A mixture with unequal weights and scales that the chain mixes well on,
  # so that the weights can be checked precisely.
  log_mix1 <- function(x) log(0.3 * stats::dnorm(x, 0, 1) + 0.7 * stats::dnorm(x, 4, 2))
  set.seed(1)
  fit <- twalk(log_mix1, n_iter = 1e5, x0 = 0, xp0 = 1, show_progress = FALSE,
               control = twalk_control(penalty_prob = 0.5))
  draws <- fit$samples[-(1:5000), 1]
  exact <- 0.3 * stats::pnorm(2, lower.tail = FALSE) +
    0.7 * stats::pnorm(2, 4, 2, lower.tail = FALSE)
  expect_equal(mean(draws > 2), exact, tolerance = 0.03)
  expect_equal(mean(draws), 2.8, tolerance = 0.05)
})

test_that("the penalised t-walk jumps between well-separated modes", {
  # Example 1 of Medina-Aguayo & Christen (2022). The standard t-walk stays in
  # the mode where it starts; the penalised one switches modes repeatedly.
  # Single chains switch slowly (the paper reports an IAT near 1000), so the
  # mode weights are checked on averages of many long chains in the vignette,
  # not here.
  log_post <- log_mixture()
  run <- function(penalty_prob) {
    set.seed(8)
    twalk(log_post, n_iter = 2e5, x0 = c(0, 0), xp0 = c(0.5, 0.5),
          show_progress = FALSE, control = twalk_control(penalty_prob = penalty_prob))
  }
  switches <- function(fit) sum(diff(fit$samples[, 1] > 10) != 0)

  expect_identical(switches(run(0)), 0L)
  expect_gt(switches(run(0.1)), 50)
})

test_that("twalk_control validates the penalised-move settings", {
  control <- twalk_control(penalty_prob = 0.1, penalty_kappa = 2)
  expect_equal(control$penalty_prob, 0.1)
  expect_equal(control$penalty_kappa, 2)
  expect_output(print(control), "penalised move: probability = 0.1, kappa = 2")
  expect_false(any(grepl("penalised", capture.output(print(twalk_control())))))

  expect_error(twalk_control(penalty_prob = 1), "`penalty_prob`")
  expect_error(twalk_control(penalty_prob = -0.1), "`penalty_prob`")
  expect_error(twalk_control(penalty_prob = c(0.1, 0.2)), "`penalty_prob`")
  expect_error(twalk_control(penalty_kappa = 0), "`penalty_kappa`")

  # Controls created by earlier versions have no penalty settings.
  legacy <- twalk_control()
  legacy$penalty_prob <- NULL
  legacy$penalty_kappa <- NULL
  expect_identical(Rtwalk:::.penalty_settings(legacy), list(prob = 0, kappa = 3))
  set.seed(9)
  fit <- twalk(log_std_normal, n_iter = 50, x0 = 0, xp0 = 1,
               show_progress = FALSE, control = legacy)
  expect_equal(nrow(fit$kernel_statistics), 4)
})
