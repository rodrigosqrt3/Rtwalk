test_that("dot_product calculates correctly", {
  expect_equal(dot_product(c(1, 2, 3), c(4, 5, 6)), 32)
})

test_that("simulate_beta covers both branch probabilities", {
  set.seed(123)
  betas <- replicate(100, simulate_beta(at = 6.0))
  expect_true(all(betas > 0))
  expect_true(length(betas) == 100)
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

test_that("twalk_move tests all 4 kernel types and both directions", {
  obj_fun <- function(p, ...) -0.5 * sum(p^2)
  supp_fun <- function(p, ...) TRUE

  x <- c(0, 0); xp <- c(1, 1)
  U <- obj_fun(x); Up <- obj_fun(xp)

  set.seed(42)
  m1 <- twalk_move(2, obj_fun, supp_fun, x, U, xp, Up, p_traverse = 1, p_walk = 0, p_blow = 0)
  expect_named(m1, c("y", "prop_U", "yp", "prop_Up", "alpha"))

  m2 <- twalk_move(2, obj_fun, supp_fun, x, U, xp, Up, p_traverse = 0, p_walk = 1, p_blow = 0)
  expect_named(m2, c("y", "prop_U", "yp", "prop_Up", "alpha"))

  m3 <- twalk_move(2, obj_fun, supp_fun, x, U, xp, Up, p_traverse = 0, p_walk = 0, p_blow = 1)
  expect_named(m3, c("y", "prop_U", "yp", "prop_Up", "alpha"))

  m4 <- twalk_move(2, obj_fun, supp_fun, x, U, xp, Up, p_traverse = 0, p_walk = 0, p_blow = 0)
  expect_named(m4, c("y", "prop_U", "yp", "prop_Up", "alpha"))

  supp_false <- function(p, ...) FALSE
  m_fail <- twalk_move(2, obj_fun, supp_false, x, U, xp, Up, p_traverse = 1, p_walk = 0, p_blow = 0)
  expect_equal(m_fail$alpha, 0)
})

test_that("twalk_move handles n_phi == 0 in Traverse kernel for both directions", {
  obj_fun <- function(p, ...) -0.5 * sum(p^2)
  supp_fun <- function(p, ...) TRUE

  x <- c(0, 0); xp <- c(1, 1)
  U <- obj_fun(x); Up <- obj_fun(xp)

  set.seed(123)
  for (i in 1:20) {
    m <- twalk_move(
      n_dim = 2,
      log_post_fun = obj_fun,
      support_fun = supp_fun,
      x = x, U = U, xp = xp, Up = Up,
      p_phi = 0,
      p_traverse = 1, p_walk = 0, p_blow = 0
    )
    expect_equal(m$alpha, 1)
  }
})
