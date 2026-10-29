# Extracted from test-kernels.R:142

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "Rtwalk", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
obj_fun <- function(p, ...) -0.5 * sum(p^2)
supp_fun <- function(p, ...) TRUE
x <- c(0, 0)
xp <- c(1, 1)
U <- obj_fun(x)
Up <- obj_fun(xp)
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
