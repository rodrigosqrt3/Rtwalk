# Extracted from test-kernels.R:97

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
set.seed(42)
m1 <- twalk_move(2, obj_fun, supp_fun, x, U, xp, Up, p_traverse = 1, p_walk = 0, p_blow = 0)
