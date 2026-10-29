test_that("single-chain results convert to coda mcmc objects", {
  log_post <- function(x) -0.5 * sum(x^2)

  set.seed(123)
  result <- twalk(
    log_post,
    n_iter = 20,
    x0 = c(-1, 1),
    xp0 = c(1, -1),
    show_progress = FALSE
  )

  primary <- coda::as.mcmc(result, burnin_frac = 0.2)
  companion <- coda::as.mcmc(
    result,
    burnin_frac = 0.2,
    companion = TRUE
  )

  expect_s3_class(primary, "mcmc")
  expect_s3_class(companion, "mcmc")
  expect_equal(
    unname(as.matrix(primary)),
    unname(result$samples[5:20, , drop = FALSE])
  )
  expect_equal(
    unname(as.matrix(companion)),
    unname(result$companion_samples[5:20, , drop = FALSE])
  )
})

test_that("parallel results convert to coda mcmc.list objects", {
  log_post <- function(x) -0.5 * sum(x^2)

  set.seed(123)
  result <- twalk(
    log_post,
    n_iter = 20,
    x0 = c(-1, 1),
    xp0 = c(1, -1),
    n_chains = 2,
    n_cores = 2,
    show_progress = FALSE
  )

  expect_error(
    coda::as.mcmc(result),
    "requires a single-chain result"
  )

  converted <- coda::as.mcmc.list(result, burnin_frac = 0.2)

  expect_s3_class(converted, "mcmc.list")
  expect_length(converted, 2)
  expect_equal(
    unname(as.matrix(converted[[1L]])),
    unname(result$individual_chains[[1L]]$samples[5:20, , drop = FALSE])
  )
  expect_equal(
    unname(as.matrix(converted[[2L]])),
    unname(result$individual_chains[[2L]]$samples[5:20, , drop = FALSE])
  )
})

test_that("coda conversion validates its options", {
  object <- structure(
    list(samples = matrix(1:10, ncol = 1)),
    class = "twalk"
  )

  expect_error(
    coda::as.mcmc(object, companion = NA),
    "`companion` must be TRUE or FALSE"
  )
})
