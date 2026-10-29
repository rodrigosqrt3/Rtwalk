test_that("simulate_beta covers both branch probabilities", {
  set.seed(123)
  betas <- replicate(100, simulate_beta(at = 6.0))
  expect_true(all(betas > 0))
  expect_true(length(betas) == 100)
})

test_that("simulate_beta follows psi_t of Christen & Fox (2010, Section 2.2)", {
  set.seed(1)
  at <- 6
  betas <- replicate(20000, simulate_beta(at))

  # P(beta <= 1) is the weight of the first mixture component, (at - 1) / (2 at)
  expect_equal(mean(betas <= 1), (at - 1) / (2 * at), tolerance = 0.02)
  # The paper reports P(beta < 2) ~ 0.98 for at = 6
  expect_equal(mean(betas < 2), 0.98, tolerance = 0.01)
  # Exact CDF at 2: 1 - (at + 1) / (2 at) * 2^(1 - at)
  expect_equal(mean(betas < 2), 1 - (at + 1) / (2 * at) * 2^(1 - at), tolerance = 0.01)
})

test_that("acceptance probabilities are stable and bounded", {
  expect_equal(acceptance_probability(log(0.25)), 0.25)
  expect_equal(acceptance_probability(0), 1)
  expect_equal(acceptance_probability(log(2)), 1)
  expect_equal(acceptance_probability(1000), 1)
  expect_equal(acceptance_probability(-Inf), 0)
  expect_equal(acceptance_probability(NaN), 0)
})

test_that("kernels and log-density helpers run without errors and cover zero-phi/small sigma branches", {
  x <- c(1, 2)
  xp <- c(1, 2)

  k_trav <- kernel_traverse(n_dim = 2, p_phi = 1.0, x = x, xp = c(2, 3), beta = 1.5)
  expect_named(k_trav, c("proposal", "n_phi"))

  k_walk <- kernel_walk(n_dim = 2, p_phi = 0.5, aw = 1.5, x = x, xp = c(2, 3))
  expect_named(k_walk, c("proposal", "n_phi"))

  k_blow <- kernel_blow(n_dim = 2, p_phi = 1.0, x = x, xp = xp)
  expect_named(k_blow, c("proposal", "n_phi", "phi"))

  expect_equal(log_density_blow(n_phi = 0, phi = c(FALSE, FALSE), h = x, x = x, xp = xp), 0)
  ld_blow <- log_density_blow(n_phi = 2, phi = c(TRUE, TRUE), h = c(1.1, 2.1), x = x, xp = xp)
  expect_true(is.numeric(ld_blow))

  k_hop <- kernel_hop(n_dim = 2, p_phi = 1.0, x = x, xp = xp)
  expect_named(k_hop, c("proposal", "n_phi", "phi"))

  expect_equal(log_density_hop(n_phi = 0, phi = c(FALSE, FALSE), h = x, x = x, xp = xp), 0)
  ld_hop <- log_density_hop(n_phi = 2, phi = c(TRUE, TRUE), h = c(1.1, 2.1), x = x, xp = xp)
  expect_true(is.numeric(ld_hop))
})

test_that("Walk proposal matches the scaled random walk of Section 2.1", {
  set.seed(5)
  aw <- 1.5
  x <- c(0.3, -1, 2)
  xp <- c(1, 0.5, -2)
  move <- kernel_walk(3, 1, aw, x, xp)

  set.seed(5)
  u <- runif(3)
  alpha <- aw / (1 + aw) * (-1 + 2 * u + aw * u^2)
  expect_equal(move$proposal, x + (x - xp) * alpha)
  expect_true(all(alpha >= -aw / (1 + aw) & alpha <= aw))
})

test_that("Blow log density only uses selected coordinates", {
  x <- 0:9
  xp <- x + 1:10
  phi <- c(rep(TRUE, 4), rep(FALSE, 6))
  n_phi <- sum(phi)
  sigma <- max(phi * abs(xp - x))

  h <- x
  h[phi] <- xp[phi] + sigma * c(0.1, -0.2, 0.3, -0.4)

  expected <-
    (n_phi / 2) * log(2 * pi) +
    n_phi * log(sigma) +
    0.5 * sum((h[phi] - xp[phi])^2) / sigma^2

  expect_equal(
    log_density_blow(n_phi, phi, h, x, xp),
    expected
  )
})

test_that("Hop log density only uses selected coordinates", {
  x <- 0:9
  xp <- x + 1:10
  phi <- c(rep(TRUE, 4), rep(FALSE, 6))
  n_phi <- sum(phi)
  sigma <- max(phi * abs(xp - x)) / 3

  h <- x
  h[phi] <- x[phi] + sigma * c(0.1, -0.2, 0.3, -0.4)
  h[!phi] <- 1000

  expected <-
    (n_phi / 2) * log(2 * pi) +
    n_phi * log(sigma) +
    0.5 * sum((h[phi] - x[phi])^2) / sigma^2

  expect_equal(
    log_density_hop(n_phi, phi, h, x, xp),
    expected
  )
})

test_that("Hop and Blow densities match the closed forms of Section 2.3", {
  x <- c(0.3, -1, 2, 5)
  xp <- c(1, 0.5, -2, 4)
  phi <- c(TRUE, FALSE, TRUE, TRUE)
  sigma <- max(abs(x - xp)[phi])

  h_hop <- x
  h_hop[phi] <- h_hop[phi] + c(0.2, -0.4, 0.1)
  g_hop <- (2 * pi)^(-3 / 2) * 3^3 / sigma^3 *
    exp(-9 / (2 * sigma^2) * sum((h_hop - x)[phi]^2))
  expect_equal(log_density_hop(3, phi, h_hop, x, xp), -log(g_hop))

  h_blow <- x
  h_blow[phi] <- xp[phi] + c(0.2, -0.4, 0.1)
  g_blow <- (2 * pi)^(-3 / 2) / sigma^3 *
    exp(-1 / (2 * sigma^2) * sum((h_blow - xp)[phi]^2))
  expect_equal(log_density_blow(3, phi, h_blow, x, xp), -log(g_blow))
})

test_that("twalk_move generates all four kernels with their Hastings terms", {
  a <- c(0, 0)
  b <- c(1, 1)

  set.seed(42)
  moves <- lapply(1:4, function(kernel) twalk_move(a, b, 1, 1.5, 6, kernel))

  for (move in moves) {
    expect_named(move, c("h", "n_phi", "log_hastings"))
    expect_length(move$h, 2)
    expect_equal(move$n_phi, 2)
  }

  expect_true(is.numeric(moves[[1]]$log_hastings))
  expect_equal(moves[[2]]$log_hastings, 0)
  expect_true(is.function(moves[[3]]$log_hastings))
  expect_true(is.function(moves[[4]]$log_hastings))
})

test_that("Traverse Jacobian is beta^(n_phi - 2)", {
  set.seed(9)
  move <- twalk_move(c(0, 0, 0), c(1, 2, 3), 1, 1.5, 6, 1L)

  set.seed(9)
  beta <- simulate_beta(6)
  expect_equal(move$log_hastings, (3 - 2) * log(beta))
  expect_equal(move$h, c(1, 2, 3) + beta * (c(1, 2, 3) - c(0, 0, 0)))
})

test_that("proposals that select no coordinate are identity moves", {
  for (kernel in 1:4) {
    set.seed(kernel)
    move <- twalk_move(c(0, 0), c(1, 1), 0, 1.5, 6, kernel)
    expect_equal(move$n_phi, 0L)
    expect_equal(move$h, c(0, 0))
  }

  log_post <- function(x) -0.5 * sum(x^2)
  set.seed(1)
  result <- twalk(
    log_post, n_iter = 200, x0 = c(0, 0), xp0 = c(1, 1),
    show_progress = FALSE, control = twalk_control(n1 = 1e-12)
  )
  expect_equal(result$acceptance_rate, 1)
  expect_equal(result$move_rate, 0)
  expect_equal(result$n_evaluations, 0)
})
