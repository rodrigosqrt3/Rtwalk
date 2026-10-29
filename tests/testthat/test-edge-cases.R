log_std_normal <- function(x) -0.5 * sum(x^2)

test_that("compiled and reference engines agree with thinning", {
  runs <- lapply(c("cpp", "R"), function(engine) {
    set.seed(31)
    twalk(log_std_normal, n_iter = 300, x0 = c(0, 0), xp0 = c(1, 1), thin = 4,
          show_progress = FALSE, control = twalk_control(engine = engine))
  })
  expect_equal(runs[[1]]$samples, runs[[2]]$samples, tolerance = 1e-10)
  expect_equal(nrow(runs[[2]]$samples), 75)
})

test_that("both engines guard their resumable interface", {
  failing <- function(x) stop("boom")
  for (engine_new in list(.engine_cpp_new, .engine_r_new)) {
    engine <- engine_new(
      x0 = c(0, 0), xp0 = c(1, 1), lp_x0 = 0, lp_xp0 = -1, n_iter = 5L,
      thin = 1L, p_phi = 1, aw = 1.5, at = 6, cum_probs = c(0, 1, 1, 1),
      store_companion = TRUE, names = NULL
    )
    expect_error(engine$reject_pending(), "no pending evaluation")
    expect_error(engine$run(failing, 1L), "boom")
    expect_true(engine$pending())
    expect_error(engine$run(log_std_normal, 1L), "pending evaluation")
    expect_equal(engine$reject_pending(), 1L)
    expect_false(engine$pending())
  }

  expect_error(twalk_engine_iter(methods::new("externalptr")), "Invalid t-walk engine")
})

test_that("invalid values returned at proposals are errors in both engines", {
  calls <- 0L
  late_vector <- function(x) {
    calls <<- calls + 1L
    if (calls > 2L) c(1, 2) else -0.5 * sum(x^2)
  }
  for (engine in c("cpp", "R")) {
    calls <- 0L
    expect_error(
      twalk(late_vector, n_iter = 10, x0 = c(0, 0), xp0 = c(1, 1),
            show_progress = FALSE, control = twalk_control(engine = engine)),
      "single numeric value"
    )
  }
})

test_that("an interruption before the first stored draw returns empty output", {
  for (engine in c("cpp", "R")) {
    calls <- 0L
    interrupting <- function(x) {
      calls <<- calls + 1L
      if (calls == 3L) {
        signalCondition(structure(
          class = c("interrupt", "condition"),
          list(message = "", call = NULL)
        ))
      }
      -0.5 * sum(x^2)
    }
    set.seed(1)
    expect_warning(
      result <- twalk(interrupting, n_iter = 50, x0 = c(0, 0), xp0 = c(1, 1),
                      show_progress = FALSE,
                      control = twalk_control(engine = engine)),
      "interrupted"
    )
    expect_equal(result$n_iter, 0L)
    expect_equal(dim(result$samples), c(0L, 2L))
    expect_equal(dim(result$companion_samples), c(0L, 2L))
  }
})

test_that("outputs without the companion trajectory stay consistent", {
  control <- twalk_control(store_companion = FALSE)
  set.seed(2)
  single <- twalk(log_std_normal, n_iter = 40, x0 = c(0, 0), xp0 = c(1, 1),
                  show_progress = FALSE, control = control)
  expect_null(single$all_samples)

  longer <- twalk_continue(single, log_std_normal, n_iter = 20, show_progress = FALSE)
  expect_equal(nrow(longer$samples), 60)
  expect_null(longer$companion_samples)

  multi <- twalk(log_std_normal, n_iter = 30,
                 x0 = rbind(c(0, 0), c(2, 2)), xp0 = rbind(c(1, 1), c(-1, -1)),
                 n_chains = 2, n_cores = 1, show_progress = FALSE, control = control)
  expect_null(multi$companion_samples)
  expect_equal(nrow(multi$samples), 60)
})

test_that("names may come from xp0 and single-row matrices are accepted", {
  set.seed(3)
  result <- twalk(log_std_normal, n_iter = 10, x0 = c(0, 0),
                  xp0 = c(a = 1, b = 1), show_progress = FALSE)
  expect_identical(colnames(result$samples), c("a", "b"))

  set.seed(3)
  from_matrix <- twalk(log_std_normal, n_iter = 10,
                       x0 = matrix(c(0, 0), 1, dimnames = list(NULL, c("a", "b"))),
                       xp0 = matrix(c(1, 1), 1), show_progress = FALSE)
  expect_equal(from_matrix$samples, result$samples)
})

test_that("chain seeds work when no random seed exists yet", {
  had_seed <- exists(".Random.seed", envir = globalenv())
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = globalenv()))

  rm(".Random.seed", envir = globalenv())
  seeds <- .chain_seeds(2)
  expect_length(seeds, 2)
  expect_false(identical(seeds[[1]], seeds[[2]]))
})

test_that("twalk_continue and twalk_step validate their inputs", {
  set.seed(4)
  fit <- twalk(log_std_normal, n_iter = 20, x0 = c(0, 0), xp0 = c(1, 1),
               thin = 5, show_progress = FALSE)

  expect_error(twalk_continue(fit, "f", 10), "must be a function")
  expect_error(twalk_continue(fit, log_std_normal, 10, n_cores = 0), "`n_cores`")
  expect_error(twalk_continue(fit, log_std_normal, 10, append = NA), "`append`")
  expect_error(twalk_continue(fit, log_std_normal, 3), "thinning interval")
  expect_error(twalk_step("f", 0, 1), "must be a function")
})

test_that("diagnostics handle constant, degenerate and very short chains", {
  constant_column <- cbind(rnorm(50), rep(3, 50))
  expect_output(table <- calculate_diagnostics(constant_column, burnin_frac = 0))
  expect_true(is.finite(table$ESS_bulk[1]))
  expect_true(all(is.na(unlist(table[2, c("ESS", "ESS_bulk", "ESS_tail", "Rhat")]))))

  short <- matrix(rnorm(6), ncol = 2)
  expect_output(table <- calculate_diagnostics(short, burnin_frac = 0))
  expect_true(all(is.na(table$ESS)))
  expect_true(all(is.na(table$Rhat)))

  # Folding around the median makes these draws constant; R-hat is undefined,
  # as in the posterior package.
  alternating <- matrix(rep(c(0, 2), 50), ncol = 1)
  expect_true(is.na(.rhat_rank(alternating)))
  skip_if_not_installed("posterior")
  expect_true(is.na(posterior::rhat(alternating)))
})

test_that("posterior conversion works for objects without names or log posterior", {
  skip_if_not_installed("posterior")
  legacy <- structure(list(samples = matrix(rnorm(40), ncol = 2)), class = "twalk")
  draws <- posterior::as_draws_array(legacy)
  expect_identical(posterior::variables(draws), c("theta1", "theta2"))
})
