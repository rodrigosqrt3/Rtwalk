# Reference implementation of the t-walk transition kernel.
#
# The functions in this file follow Christen & Fox (2010) as literally as
# possible and are the specification for the compiled engine in
# src/twalk_engine.cpp. Both engines draw random numbers in the same order, so
# they produce identical chains for the same seed (see tests/testthat).

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

# Convert a log Metropolis--Hastings ratio to the acceptance probability
# min(1, exp(log_ratio)) without overflowing for large positive ratios.
acceptance_probability <- function(log_ratio) {
  if (!is.numeric(log_ratio) || length(log_ratio) != 1L || is.na(log_ratio)) {
    return(0.0)
  }

  if (log_ratio >= 0) {
    return(1.0)
  }

  exp(log_ratio)
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
    # The continuous part of the Blow proposal is defined only on the
    # coordinates selected by phi. The remaining coordinates are point masses
    # at x and must not contribute to the Gaussian density.
    diff_vec <- h[phi] - xp[phi]
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
    # The continuous part of the Hop proposal is defined only on the
    # coordinates selected by phi. The remaining coordinates are point masses
    # at x and must not contribute to the Gaussian density.
    diff_vec <- h[phi] - x[phi]
    (n_phi / 2) * log(2 * pi) + n_phi * log(sigma) +
      0.5 * sum(diff_vec^2) / (sigma^2)
  } else {
    0
  }
}

#' Propose a single t-walk move
#'
#' Selects one of the four kernels and one of the two points, generates a
#' proposal and returns everything needed to compute the Metropolis--Hastings
#' acceptance probability. Following equation (1) of Christen & Fox (2010),
#' only one point moves: `a` is the point being moved and `b` the other one,
#' so the proposal is `h(a, b)`.
#'
#' @param a Point being moved.
#' @param b The other point.
#' @param p_phi Probability of updating each coordinate.
#' @param aw,at Walk and Traverse parameters.
#' @param kernel Integer kernel index (1 = traverse, 2 = walk, 3 = blow,
#'   4 = hop).
#' @return A list with the proposal `h`, `n_phi`, and the log Hastings
#'   correction `log_hastings` (excluding the target ratio).
#' @keywords internal
twalk_move <- function(a, b, p_phi, aw, at, kernel) {
  n_dim <- length(a)

  if (kernel == 1L) {
    beta <- simulate_beta(at)
    res <- kernel_traverse(n_dim, p_phi, a, b, beta)
    return(list(
      h = res$proposal,
      n_phi = res$n_phi,
      log_hastings = (res$n_phi - 2) * log(beta)
    ))
  }

  if (kernel == 2L) {
    res <- kernel_walk(n_dim, p_phi, aw, a, b)
    return(list(h = res$proposal, n_phi = res$n_phi, log_hastings = 0))
  }

  if (kernel == 3L) {
    res <- kernel_blow(n_dim, p_phi, a, b)
    h <- res$proposal
    hastings <- function() {
      log_density_blow(res$n_phi, res$phi, h, a, b) -
        log_density_blow(res$n_phi, res$phi, a, h, b)
    }
  } else {
    res <- kernel_hop(n_dim, p_phi, a, b)
    h <- res$proposal
    hastings <- function() {
      log_density_hop(res$n_phi, res$phi, h, a, b) -
        log_density_hop(res$n_phi, res$phi, a, h, b)
    }
  }

  list(h = h, n_phi = res$n_phi, log_hastings = hastings)
}

.kernel_names <- c("traverse", "walk", "blow", "hop")

.check_log_posterior_value <- function(value) {
  if (!is.numeric(value) || length(value) != 1L) {
    stop("`log_posterior` must return a single numeric value.", call. = FALSE)
  }

  value <- as.numeric(value)
  if (is.na(value) || (is.infinite(value) && value > 0)) {
    stop(
      "`log_posterior` must return a finite value or -Inf outside the support.",
      call. = FALSE
    )
  }

  value
}

# Reference engine. It exposes the same resumable interface as the compiled
# engine: `run()` may be interrupted by an error raised by the log posterior,
# after which `pending()` is TRUE and `reject_pending()` completes the
# iteration as a rejection.
.engine_r_new <- function(x0, xp0, lp_x0, lp_xp0, n_iter, thin, p_phi,
                          aw, at, cum_probs, store_companion, names) {
  state <- new.env(parent = emptyenv())
  n_dim <- length(x0)
  n_store <- n_iter %/% thin

  state$x <- as.numeric(x0)
  state$xp <- as.numeric(xp0)
  state$lp_x <- lp_x0
  state$lp_xp <- lp_xp0
  state$iter <- 0L
  state$n_evals <- 0L
  state$pending <- FALSE
  state$out_x <- matrix(0, n_store, n_dim)
  state$out_lp_x <- numeric(n_store)
  if (store_companion) {
    state$out_xp <- matrix(0, n_store, n_dim)
    state$out_lp_xp <- numeric(n_store)
  }
  state$proposed <- integer(4)
  state$accepted <- integer(4)
  state$moved <- integer(4)

  store_iteration <- function() {
    state$iter <- state$iter + 1L
    if (state$iter %% thin != 0L) {
      return(invisible(NULL))
    }
    row <- state$iter %/% thin
    state$out_x[row, ] <- state$x
    state$out_lp_x[row] <- state$lp_x
    if (store_companion) {
      state$out_xp[row, ] <- state$xp
      state$out_lp_xp[row] <- state$lp_xp
    }
    invisible(NULL)
  }

  evaluate <- function(log_posterior, h) {
    params <- h
    names(params) <- names
    state$pending <- TRUE
    value <- log_posterior(params)
    state$pending <- FALSE
    state$n_evals <- state$n_evals + 1L
    .check_log_posterior_value(value)
  }

  step <- function(log_posterior) {
    u_kernel <- stats::runif(1)
    kernel <- match(TRUE, u_kernel < cum_probs[1:3], nomatch = 4L)
    move_xp <- stats::runif(1) < 0.5

    if (move_xp) {
      a <- state$xp
      lp_a <- state$lp_xp
      b <- state$x
    } else {
      a <- state$x
      lp_a <- state$lp_x
      b <- state$xp
    }

    move <- twalk_move(a, b, p_phi, aw, at, kernel)
    h <- move$h
    state$proposed[kernel] <- state$proposed[kernel] + 1L

    alpha <- 0
    lp_h <- lp_a
    if (move$n_phi == 0L) {
      alpha <- 1
    } else if (kernel == 1L || any(h != b)) {
      lp_h <- evaluate(log_posterior, h)
      if (is.finite(lp_h)) {
        log_hastings <- move$log_hastings
        if (is.function(log_hastings)) log_hastings <- log_hastings()
        alpha <- acceptance_probability((lp_h - lp_a) + log_hastings)
      }
    }

    if (stats::runif(1) < alpha) {
      changed <- any(h != a)
      if (move_xp) {
        state$xp <- h
        state$lp_xp <- lp_h
      } else {
        state$x <- h
        state$lp_x <- lp_h
      }
      state$accepted[kernel] <- state$accepted[kernel] + 1L
      if (changed) state$moved[kernel] <- state$moved[kernel] + 1L
    }

    store_iteration()
  }

  list(
    run = function(log_posterior, n_steps) {
      if (state$pending) {
        stop("The t-walk engine has a pending evaluation; reject it first.", call. = FALSE)
      }
      last <- min(n_iter, state$iter + n_steps)
      while (state$iter < last) step(log_posterior)
      state$iter
    },
    reject_pending = function() {
      if (!state$pending) {
        stop("The t-walk engine has no pending evaluation.", call. = FALSE)
      }
      state$pending <- FALSE
      stats::runif(1)
      store_iteration()
      state$iter
    },
    pending = function() state$pending,
    iter = function() state$iter,
    result = function() {
      n_done <- state$iter %/% thin
      keep <- seq_len(n_done)
      list(
        samples = state$out_x[keep, , drop = FALSE],
        companion_samples = if (store_companion) state$out_xp[keep, , drop = FALSE],
        log_posterior = state$out_lp_x[keep],
        companion_log_posterior = if (store_companion) state$out_lp_xp[keep],
        proposed = state$proposed,
        accepted = state$accepted,
        moved = state$moved,
        x = state$x,
        xp = state$xp,
        lp_x = state$lp_x,
        lp_xp = state$lp_xp,
        iter = state$iter,
        n_evals = state$n_evals
      )
    }
  )
}

.engine_cpp_new <- function(x0, xp0, lp_x0, lp_xp0, n_iter, thin, p_phi,
                            aw, at, cum_probs, store_companion, names) {
  ptr <- twalk_engine_new(
    as.numeric(x0), as.numeric(xp0), lp_x0, lp_xp0,
    as.integer(n_iter), as.integer(thin), p_phi, aw, at,
    as.numeric(cum_probs), store_companion, names
  )

  list(
    run = function(log_posterior, n_steps) {
      twalk_engine_run(ptr, log_posterior, as.integer(n_steps))
    },
    reject_pending = function() twalk_engine_reject_pending(ptr),
    pending = function() twalk_engine_pending(ptr),
    iter = function() twalk_engine_iter(ptr),
    result = function() twalk_engine_result(ptr)
  )
}
