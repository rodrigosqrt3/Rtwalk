# Basic helper functions
dot_product <- function(x, y) { sum(x * y) }

#' Simulate the beta parameter for the Traverse kernel
#'
#' Samples beta from the proposal distribution in Section 2.2 of
#' Christen & Fox (2010).
#' @param at Shape parameter `a_t` (default 6).
#' @return A numeric value for beta.
#' @keywords internal
simulate_beta <- function(at = 6.0) {
  if (stats::runif(1) < (at - 1) / (2 * at)) {
    exp(1 / (at + 1) * log(stats::runif(1)))
  } else {
    exp(1 / (1 - at) * log(stats::runif(1)))
  }
}

#' Generate a proposal for the Traverse kernel
#' @param n_dim Dimension of the parameter space.
#' @param p_phi Probability of updating each coordinate.
#' @param x Current point `x`.
#' @param xp Current point `x'`.
#' @param beta Step parameter beta.
#' @return A list with the `proposal` and `n_phi` (number of moved coordinates).
#' @keywords internal
kernel_traverse <- function(n_dim, p_phi, x, xp, beta) {
  phi <- (stats::runif(n_dim) < p_phi)
  proposal <- xp + beta * (xp - x)
  proposal[!phi] <- x[!phi]
  list(proposal = proposal, n_phi = sum(phi))
}

#' Generate a proposal for the Walk kernel
#' @param n_dim Dimension of the parameter space.
#' @param p_phi Probability of updating each coordinate.
#' @param aw Scale parameter `a_w` (default 1.5).
#' @param x Current point `x`.
#' @param xp Current point `x'`.
#' @return A list with the `proposal` and `n_phi`.
#' @keywords internal
kernel_walk <- function(n_dim, p_phi, aw, x, xp) {
  u <- stats::runif(n_dim)
  phi <- (stats::runif(n_dim) < p_phi)
  z <- (aw / (1 + aw)) * (aw * u^2 + 2 * u - 1)
  z <- z * phi
  list(proposal = x + (x - xp) * z, n_phi = sum(phi))
}

#' Generate a proposal for the Blow kernel
#' @param n_dim Dimension of the parameter space.
#' @param p_phi Probability of updating each coordinate.
#' @param x Current point `x`.
#' @param xp Current point `x'`.
#' @return A list with the `proposal`, `n_phi`, and `phi`.
#' @keywords internal
kernel_blow <- function(n_dim, p_phi, x, xp) {
  phi <- (stats::runif(n_dim) < p_phi)
  sigma <- max(phi * abs(xp - x))
  if (sigma < .Machine$double.eps) sigma <- .Machine$double.eps # Avoid sigma = 0

  noise <- stats::rnorm(n_dim) * phi
  proposal <- xp * phi + sigma * noise + x * (1 - phi)

  list(proposal = proposal, n_phi = sum(phi), phi = phi)
}

#' Calculate the log of the proposal density for the Blow kernel
#' @param n_phi Number of moved coordinates.
#' @param phi Boolean vector indicating moved coordinates.
#' @param h Proposed point.
#' @param x Reference point.
#' @param xp Reference point.
#' @return The value of -log(g(h)).
#' @keywords internal
log_density_blow <- function(n_phi, phi, h, x, xp) {
  sigma <- max(phi * abs(xp - x))
  if (sigma < .Machine$double.eps) sigma <- .Machine$double.eps

  if (n_phi > 0) {
    diff_vec <- h - xp
    (n_phi / 2) * log(2 * pi) + n_phi * log(sigma) + 0.5 * sum(diff_vec^2) / (sigma^2)
  } else {
    0
  }
}

#' Generate a proposal for the Hop kernel
#' @param n_dim Dimension of the parameter space.
#' @param p_phi Probability of updating each coordinate.
#' @param x Current point `x`.
#' @param xp Current point `x'`.
#' @return A list with the `proposal`, `n_phi`, and `phi`.
#' @keywords internal
kernel_hop <- function(n_dim, p_phi, x, xp) {
  phi <- (stats::runif(n_dim) < p_phi)
  sigma <- max(phi * abs(xp - x)) / 3
  if (sigma < .Machine$double.eps) sigma <- .Machine$double.eps

  noise <- stats::rnorm(n_dim) * phi
  proposal <- x + sigma * noise
  proposal[!phi] <- x[!phi]

  list(proposal = proposal, n_phi = sum(phi), phi = phi)
}

#' Calculate the log of the proposal density for the Hop kernel
#' @param n_phi Number of moved coordinates.
#' @param phi Boolean vector indicating moved coordinates.
#' @param h Proposed point.
#' @param x Reference point.
#' @param xp Reference point.
#' @return The value of -log(g(h)).
#' @keywords internal
log_density_hop <- function(n_phi, phi, h, x, xp) {
  sigma <- max(phi * abs(xp - x)) / 3
  if (sigma < .Machine$double.eps) sigma <- .Machine$double.eps

  if (n_phi > 0) {
    diff_vec <- h - x
    (n_phi / 2) * log(2*pi) - n_phi * log(3) + n_phi * log(max(phi * abs(xp - x))) + 0.5 * sum(diff_vec^2) / (sigma^2)
  } else {
    0
  }
}

#' Run a single t-walk move (step)
#'
#' Selects one of the four kernels, generates a proposal, and calculates
#' the Metropolis-Hastings acceptance probability.
#' @param ... Arguments passed to `log_post_fun` and `support_fun`.
#' @return A list containing the proposal and the acceptance probability.
#' @keywords internal
twalk_move <- function(n_dim, log_post_fun, support_fun, x, U, xp, Up,
                           at = 6.0, aw = 1.5,
                           p_phi = min(n_dim, 4) / n_dim,
                           p_traverse = 0.4918, p_walk = 0.4918, p_blow = 0.0082, ...) {

  # Kernel selection probabilities
  p_hop <- 1 - p_traverse - p_walk - p_blow
  kernel_probs <- c(p_traverse, p_walk, p_blow, p_hop)

  kernel_choice <- sample.int(4, 1, prob = kernel_probs)
  direction <- stats::runif(1)

  # Initialize proposal variables
  y <- NULL; yp <- NULL; prop_U <- NULL; prop_Up <- NULL; alpha <- 0.0

  if (kernel_choice == 1) { # Traverse Kernel
    beta <- simulate_beta(at)
    if (direction < 0.5) {
      res <- kernel_traverse(n_dim, p_phi, xp, x, beta)
      yp <- res$proposal; n_phi <- res$n_phi; y <- x; prop_U <- U
      if (support_fun(yp, ...)) {
        prop_Up <- log_post_fun(yp, ...); if (n_phi == 0) alpha <- 1 else alpha <- exp((U - prop_U) + (Up - prop_Up) + (n_phi - 2) * log(beta))
      }
    } else {
      res <- kernel_traverse(n_dim, p_phi, x, xp, beta)
      y <- res$proposal; n_phi <- res$n_phi; yp <- xp; prop_Up <- Up
      if (support_fun(y, ...)) {
        prop_U <- log_post_fun(y, ...); if (n_phi == 0) alpha <- 1 else alpha <- exp((U - prop_U) + (Up - prop_Up) + (n_phi - 2) * log(beta))
      }
    }
  } else if (kernel_choice == 2) { # Walk Kernel
    if (direction < 0.5) {
      res <- kernel_walk(n_dim, p_phi, aw, xp, x)
      yp <- res$proposal; y <- x; prop_U <- U
      if (support_fun(yp, ...) && any(yp != y)) {
        prop_Up <- log_post_fun(yp, ...); alpha <- exp((U - prop_U) + (Up - prop_Up))
      }
    } else {
      res <- kernel_walk(n_dim, p_phi, aw, x, xp)
      y <- res$proposal; yp <- xp; prop_Up <- Up
      if (support_fun(y, ...) && any(y != yp)) {
        prop_U <- log_post_fun(y, ...); alpha <- exp((U - prop_U) + (Up - prop_Up))
      }
    }
  } else if (kernel_choice == 3) { # Blow Kernel
    if (direction < 0.5) {
      res <- kernel_blow(n_dim, p_phi, xp, x)
      yp <- res$proposal; n_phi <- res$n_phi; phi <- res$phi; y <- x; prop_U <- U
      if (support_fun(yp, ...) && any(yp != x)) {
        prop_Up <- log_post_fun(yp, ...); W1 <- log_density_blow(n_phi, phi, yp, xp, x); W2 <- log_density_blow(n_phi, phi, xp, yp, x); alpha <- exp((U - prop_U) + (Up - prop_Up) + (W1 - W2))
      }
    } else {
      res <- kernel_blow(n_dim, p_phi, x, xp)
      y <- res$proposal; n_phi <- res$n_phi; phi <- res$phi; yp <- xp; prop_Up <- Up
      if (support_fun(y, ...) && any(y != xp)) {
        prop_U <- log_post_fun(y, ...); W1 <- log_density_blow(n_phi, phi, y, x, xp); W2 <- log_density_blow(n_phi, phi, x, y, xp); alpha <- exp((U - prop_U) + (Up - prop_Up) + (W1 - W2))
      }
    }
  } else { # Hop Kernel
    if (direction < 0.5) {
      res <- kernel_hop(n_dim, p_phi, xp, x)
      yp <- res$proposal; n_phi <- res$n_phi; phi <- res$phi; y <- x; prop_U <- U
      if (support_fun(yp, ...) && any(yp != x)) {
        prop_Up <- log_post_fun(yp, ...); W1 <- log_density_hop(n_phi, phi, yp, xp, x); W2 <- log_density_hop(n_phi, phi, xp, yp, x); alpha <- exp((U - prop_U) + (Up - prop_Up) + (W1 - W2))
      }
    } else {
      res <- kernel_hop(n_dim, p_phi, x, xp)
      y <- res$proposal; n_phi <- res$n_phi; phi <- res$phi; yp <- xp; prop_Up <- Up
      if (support_fun(y, ...) && any(y != xp)) {
        prop_U <- log_post_fun(y, ...); W1 <- log_density_hop(n_phi, phi, y, x, xp); W2 <- log_density_hop(n_phi, phi, x, y, xp); alpha <- exp((U - prop_U) + (Up - prop_Up) + (W1 - W2))
      }
    }
  }

  if (is.nan(alpha)) { alpha <- 0.0 }

  list(y = y, prop_U = prop_U, yp = yp, prop_Up = prop_Up, alpha = alpha)
}
