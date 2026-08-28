.detect_available_cores <- function(detector = parallel::detectCores) {
  detected_cores <- detector()
  if (is.na(detected_cores)) {
    detected_cores <- 1L
  }

  detected_cores
}

#' Run the t-walk MCMC Algorithm
#'
#' This function implements the t-walk algorithm by Christen & Fox (2010),
#' a general-purpose MCMC sampler that does not require manual tuning.
#' The function can run multiple independent MCMC chains in parallel
#' to accelerate execution and facilitate convergence diagnostics.
#'
#' @param log_posterior A function that takes a parameter vector as its
#'   first argument and returns one numeric log posterior density. It may
#'   return `-Inf` outside the support, but must not return vectors, `NA`,
#'   `NaN`, or `+Inf`.
#'   Additional arguments can be passed to this function via `...`.
#' @param n_iter The number of iterations to run for each chain.
#' @param x0 A numeric vector with the initial values for the first point (`x`).
#' @param xp0 A numeric vector with the initial values for the second point (`x'`).
#' @param n_chains The number of independent MCMC chains to run.
#'   Defaults to `1`, which runs a single chain sequentially. If greater
#'   than 1, parallel mode is activated.
#' @param n_cores The number of CPU cores to use in parallel mode.
#'   If `NULL` (default), it will attempt to use all available cores minus one.
#'   Parallel random-number streams are initialized from R's current random
#'   state, so calling `set.seed()` before `twalk()` makes results reproducible.
#' @param show_progress Logical; whether to display progress bars and status
#'   messages. Defaults to `TRUE`.
#' @param ... Additional arguments to be passed to the `log_posterior` function.
#'
#' @return A list containing:
#' \item{samples}{The primary t-walk trajectory, with `n_iter` rows per chain.}
#' \item{companion_samples}{The auxiliary trajectory maintained by the t-walk.}
#' \item{all_samples}{Legacy concatenation of the primary and auxiliary
#'       trajectories. This is retained for compatibility and should not be
#'       treated as a single time-ordered MCMC chain.}
#' \item{acceptance_rate}{The average Metropolis--Hastings acceptance rate
#'       across all chains. Accepted identity proposals are included, as
#'       required by the t-walk transition kernel.}
#' \item{move_rate}{The proportion of iterations in which an accepted
#'       proposal actually changed at least one of the two t-walk points.}
#' \item{no_move_rate}{The proportion of iterations containing an accepted
#'       identity proposal. It equals `acceptance_rate - move_rate`.}
#' \item{n_iter}{The number of iterations generated per chain.}
#' \item{n_chains}{The number of independent chains.}
#' \item{total_iterations}{The total number of primary samples generated
#'       (`n_iter * n_chains`).}
#' \item{n_dim}{The dimension of the parameter space.}
#' \item{individual_chains}{If `n_chains > 1`, a list containing the raw
#'       results from each separate chain, useful for diagnostics like R-hat.}
#'
#' @export
#' @importFrom parallel detectCores makeCluster clusterEvalQ clusterExport clusterSetRNGStream parLapply stopCluster
#' @importFrom stats rnorm runif
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @examples
#' # Example 1: Sampling from a Bivariate Normal (sequential mode)
#' # The 'mvtnorm' package is required for this example
#' if (requireNamespace("mvtnorm", quietly = TRUE)) {
#'   log_post <- function(x) {
#'     mvtnorm::dmvnorm(x, mean = c(0, 0), sigma = matrix(c(1, 0.8, 0.8, 1), 2, 2), log = TRUE)
#'   }
#'
#'   # Run with fewer iterations for a quick example
#'   # Set a seed for reproducibility
#'   set.seed(123)
#'   result_seq <- twalk(log_posterior = log_post, n_iter = 5000,
#'                           x0 = c(-1, 1), xp0 = c(1, -1))
#'
#'   plot(result_seq$samples, pch = '.', main = "t-walk Samples (Sequential)")
#' }
#'
#' \donttest{
#' # Example 2: The same problem in parallel (will run faster)
#' # Using 2 chains. n_iter is now per chain.
#' if (requireNamespace("mvtnorm", quietly = TRUE)) {
#'   set.seed(123)
#'   result_par <- twalk(log_posterior = log_post, n_iter = 2500,
#'                           x0 = c(-1, 1), xp0 = c(1, -1), n_chains = 2)
#'
#'   plot(result_par$samples, pch = '.', main = "t-walk Samples (Parallel)")
#' }
#' }
twalk <- function(log_posterior, n_iter, x0, xp0,
                  n_chains = 1, n_cores = NULL, show_progress = TRUE, ...) {

  is_positive_integer <- function(x) {
    is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
      x >= 1 && x == floor(x) && x <= .Machine$integer.max
  }

  if (!is.function(log_posterior)) {
    stop("`log_posterior` must be a function.", call. = FALSE)
  }

  if (!is_positive_integer(n_iter)) {
    stop("`n_iter` must be a positive integer.", call. = FALSE)
  }

  if (!is_positive_integer(n_chains)) {
    stop("`n_chains` must be a positive integer.", call. = FALSE)
  }

  if (!is.null(n_cores) && !is_positive_integer(n_cores)) {
    stop("`n_cores` must be NULL or a positive integer.", call. = FALSE)
  }

  if (!is.logical(show_progress) || length(show_progress) != 1L || is.na(show_progress)) {
    stop("`show_progress` must be TRUE or FALSE.", call. = FALSE)
  }

  valid_initial_point <- function(x) {
    is.numeric(x) && is.null(dim(x)) && length(x) > 0L && all(is.finite(x))
  }

  if (!valid_initial_point(x0) || !valid_initial_point(xp0)) {
    stop("`x0` and `xp0` must be non-empty numeric vectors with finite values.", call. = FALSE)
  }

  if (length(x0) != length(xp0)) {
    stop("`x0` and `xp0` must have the same length.", call. = FALSE)
  }

  n_iter <- as.integer(n_iter)
  n_chains <- as.integer(n_chains)
  if (!is.null(n_cores)) {
    n_cores <- as.integer(n_cores)
  }

  # Capture all extra arguments in a list
  extra_args <- list(...)

  # --- SEQUENTIAL BLOCK ---
  if (n_chains == 1) {
    if (show_progress) {
      message("--- Running t-walk in sequential mode (1 chain) ---")
    }

    n_dim <- length(x0)

    internal_args <- extra_args

    evaluation_cache <- new.env(parent = emptyenv())
    evaluation_cache$has_value <- FALSE

    evaluate_log_posterior <- function(params) {
      if (evaluation_cache$has_value &&
          identical(params, evaluation_cache$params)) {
        return(evaluation_cache$value)
      }

      res <- tryCatch(
        do.call(log_posterior, c(list(params), internal_args)),
        error = function(e) -Inf
      )

      if (!is.numeric(res) || length(res) != 1L) {
        stop("`log_posterior` must return a single numeric value.", call. = FALSE)
      }

      if (is.na(res) || is.nan(res) || (is.infinite(res) && res > 0)) {
        stop("`log_posterior` must return a finite value or -Inf outside the support.", call. = FALSE)
      }

      evaluation_cache$params <- params
      evaluation_cache$value <- res
      evaluation_cache$has_value <- TRUE

      res
    }

    # Wrapper for the objective function (-log_posterior)
    objective_fun <- function(params, ...) {
      -evaluate_log_posterior(params)
    }

    # Wrapper for the support function
    support_fun <- function(params, ...) {
      is.finite(evaluate_log_posterior(params))
    }

    initial_log_x <- evaluate_log_posterior(x0)
    initial_log_xp <- evaluate_log_posterior(xp0)

    if (!is.finite(initial_log_x) || !is.finite(initial_log_xp)) {
      stop("Initial points are outside the support (log-posterior is -Inf or returns an error).")
    }

    U <- -initial_log_x
    Up <- -initial_log_xp
    x_current <- x0; xp_current <- xp0

    x_samples <- matrix(NA, nrow = n_iter, ncol = n_dim)
    xp_samples <- matrix(NA, nrow = n_iter, ncol = n_dim)
    accepted_count <- 0L
    moved_count <- 0L

    use_progress_bar <- show_progress
    if (use_progress_bar) {
      progress_bar <- utils::txtProgressBar(min = 0, max = n_iter, style = 3, width = 50, char = "=")
    }

    for (i in 1:n_iter) {
      move <- twalk_move(
        n_dim = n_dim,
        log_post_fun = objective_fun,
        support_fun = support_fun,
        x = x_current,
        U = U,
        xp = xp_current,
        Up = Up
      )

      if (stats::runif(1) < move$alpha) {
        state_changed <- any(move$y != x_current) ||
          any(move$yp != xp_current)

        x_current <- move$y
        U <- move$prop_U
        xp_current <- move$yp
        Up <- move$prop_Up
        accepted_count <- accepted_count + 1L
        if (state_changed) {
          moved_count <- moved_count + 1L
        }
      }

      x_samples[i, ] <- x_current
      xp_samples[i, ] <- xp_current
      if (use_progress_bar) {
        utils::setTxtProgressBar(progress_bar, i)
      }
    }
    if (use_progress_bar) {
      close(progress_bar)
    }

    acceptance_rate <- accepted_count / n_iter
    move_rate <- moved_count / n_iter
    no_move_rate <- (accepted_count - moved_count) / n_iter
    if (use_progress_bar) {
      message(sprintf("\nAcceptance rate: %.2f%%", acceptance_rate * 100))
      message(sprintf("Actual move rate: %.2f%%", move_rate * 100))
    }

    return(structure(
      list(
        samples = x_samples,
        companion_samples = xp_samples,
        all_samples = rbind(x_samples, xp_samples),
        acceptance_rate = acceptance_rate,
        move_rate = move_rate,
        no_move_rate = no_move_rate,
        n_iter = n_iter,
        n_chains = 1L,
        total_iterations = n_iter,
        n_dim = n_dim
      ),
      class = "twalk"
    ))
  }

  # --- PARALLEL BLOCK ---
  else {

    if (is.null(n_cores)) {
      detected_cores <- .detect_available_cores()
      n_cores <- max(1L, detected_cores - 1L)
    }
    n_cores_used <- min(n_chains, n_cores)

    if (show_progress) {
      message(sprintf("--- Running t-walk in PARALLEL mode (%d chains on %d cores) ---", n_chains, n_cores_used))
    }

    cl <- parallel::makeCluster(n_cores_used)
    on.exit(parallel::stopCluster(cl))

    # Derive reproducible, independent worker streams from the RNG state in the
    # main R process. This respects set.seed() without resetting the user's RNG.
    parallel_seed <- sample.int(.Machine$integer.max, 1)
    parallel::clusterSetRNGStream(cl, iseed = parallel_seed)

    # Export the log_posterior function and extra_args to cluster workers
    parallel::clusterExport(cl, varlist = c("log_posterior", "extra_args", "x0", "xp0", "n_iter"), envir = environment())

    # Functions created in the global or a local environment may refer to data
    # and helper functions stored alongside them. PSOCK workers start with a
    # clean workspace, so recursively collect and export those dependencies.
    collect_dependencies <- function(fun, collected = list()) {
      fun_env <- environment(fun)
      if (is.null(fun_env)) {
        return(collected)
      }

      global_names <- codetools::findGlobals(fun, merge = TRUE)
      local_names <- global_names[vapply(
        global_names,
        exists,
        logical(1),
        envir = fun_env,
        inherits = FALSE
      )]

      for (name in setdiff(local_names, names(collected))) {
        object <- get(name, envir = fun_env, inherits = FALSE)
        collected[[name]] <- object

        if (is.function(object)) {
          collected <- collect_dependencies(object, collected)
        }
      }

      collected
    }

    dependencies <- collect_dependencies(log_posterior)
    if (length(dependencies) > 0L) {
      dependency_env <- list2env(dependencies, parent = emptyenv())
      parallel::clusterExport(
        cl,
        varlist = names(dependencies),
        envir = dependency_env
      )
    }

    # Load required packages on each worker node
    parallel::clusterEvalQ(cl, {
      # Add any packages your log_posterior might need
      # e.g., library(mvtnorm)
    })

    # This is the function that will be executed on each worker node
    run_single_chain <- function(chain_index) {
      # Use 'do.call' to safely construct the function call,
      # passing the extra arguments (...) correctly.
      chain_result <- do.call(twalk, c(
        list(log_posterior = log_posterior, n_iter = n_iter, x0 = x0, xp0 = xp0,
             n_chains = 1, show_progress = FALSE),
        extra_args
      ))
      return(chain_result)
    }

    if (show_progress) {
      message("Distributing work among cores...")
    }
    results_list <- parallel::parLapply(cl, 1:n_chains, run_single_chain)

    if (show_progress) {
      message("Chains completed. Combining results...")
    }

    combined_primary_samples <- do.call(rbind, lapply(results_list, function(res) res$samples))
    combined_companion_samples <- do.call(rbind, lapply(results_list, function(res) res$companion_samples))
    combined_samples <- do.call(rbind, lapply(results_list, function(res) res$all_samples))
    chain_acceptance_rates <- vapply(
      results_list,
      function(res) res$acceptance_rate,
      numeric(1)
    )
    mean_acceptance_rate <- mean(chain_acceptance_rates)

    # Derive actual movements from the returned trajectories instead of
    # relying on additional fields in the worker result. Besides being an
    # independent consistency check, this also keeps PSOCK execution robust
    # when a development version is tested against an older installed build.
    chain_move_rates <- vapply(results_list, function(res) {
      previous_primary <- rbind(
        x0,
        res$samples[-nrow(res$samples), , drop = FALSE]
      )
      previous_companion <- rbind(
        xp0,
        res$companion_samples[-nrow(res$companion_samples), , drop = FALSE]
      )

      primary_changed <- rowSums(res$samples != previous_primary) > 0L
      companion_changed <- rowSums(
        res$companion_samples != previous_companion
      ) > 0L

      mean(primary_changed | companion_changed)
    }, numeric(1))

    for (chain_index in seq_along(results_list)) {
      results_list[[chain_index]]$move_rate <- chain_move_rates[[chain_index]]
      results_list[[chain_index]]$no_move_rate <-
        chain_acceptance_rates[[chain_index]] - chain_move_rates[[chain_index]]
    }

    mean_move_rate <- mean(chain_move_rates)
    mean_no_move_rate <- mean_acceptance_rate - mean_move_rate
    if (show_progress) {
      message(sprintf("\nMean acceptance rate across chains: %.2f%%", mean_acceptance_rate * 100))
      message(sprintf("Mean actual move rate across chains: %.2f%%", mean_move_rate * 100))
    }

    return(structure(
      list(
        samples = combined_primary_samples,
        companion_samples = combined_companion_samples,
        all_samples = combined_samples,
        acceptance_rate = mean_acceptance_rate,
        move_rate = mean_move_rate,
        no_move_rate = mean_no_move_rate,
        n_iter = n_iter,
        n_chains = n_chains,
        total_iterations = n_iter * n_chains,
        n_dim = length(x0),
        individual_chains = results_list
      ),
      class = "twalk"
    ))
  }
}
