#' Run the t-walk MCMC Algorithm
#'
#' This function implements the t-walk algorithm by Christen & Fox (2010),
#' a general-purpose MCMC sampler that does not require manual tuning.
#' The function can run multiple independent MCMC chains in parallel
#' to accelerate execution and facilitate convergence diagnostics.
#'
#' @param log_posterior A function that takes a parameter vector as its
#'   first argument and returns the scalar log posterior density.
#'   Additional arguments can be passed to this function via `...`.
#' @param n_iter The number of iterations to run for each chain.
#' @param x0 A numeric vector with the initial values for the first point (`x`).
#' @param xp0 A numeric vector with the initial values for the second point (`x'`).
#' @param n_chains The number of independent MCMC chains to run.
#'   Defaults to `1`, which runs a single chain sequentially. If greater
#'   than 1, parallel mode is activated.
#' @param n_cores The number of CPU cores to use in parallel mode.
#'   If `NULL` (default), it will attempt to use all available cores minus one.
#' @param ... Additional arguments to be passed to the `log_posterior` function.
#'
#' @return A list containing:
#' \item{all_samples}{A matrix with the combined samples from all chains.}
#' \item{acceptance_rate}{The average acceptance rate across all chains.}
#' \item{total_iterations}{The total number of samples generated (n_iter * n_chains).}
#' \item{n_dim}{The dimension of the parameter space.}
#' \item{individual_chains}{If `n_chains > 1`, a list containing the raw
#'       results from each separate chain, useful for diagnostics like R-hat.}
#'
#' @export
#' @importFrom parallel detectCores makeCluster clusterEvalQ clusterExport parLapply stopCluster
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
#'   plot(result_seq$all_samples, pch = '.', main = "t-walk Samples (Sequential)")
#' }
#'
#' \dontrun{
#' # Example 2: The same problem in parallel (will run faster)
#' # Using 2 chains. n_iter is now per chain.
#' if (requireNamespace("mvtnorm", quietly = TRUE)) {
#'   set.seed(123)
#'   result_par <- twalk(log_posterior = log_post, n_iter = 2500,
#'                           x0 = c(-1, 1), xp0 = c(1, -1), n_chains = 2)
#'
#'   plot(result_par$all_samples, pch = '.', main = "t-walk Samples (Parallel)")
#' }
#' }
twalk <- function(log_posterior, n_iter, x0, xp0,
                      n_chains = 1, n_cores = NULL, ...) {

  # Capture all extra arguments in a list
  extra_args <- list(...)

  # --- SEQUENTIAL BLOCK ---
  if (n_chains == 1) {

    is_internal_call <- "internal_call" %in% names(extra_args)

    if (!is_internal_call) {
      cat("--- Running t-walk in sequential mode (1 chain) ---\n")
    }

    n_dim <- length(x0)

    # Create a clean copy of extra arguments for internal use,
    # removing the 'internal_call' flag.
    internal_args <- extra_args
    if (is_internal_call) {
      internal_args$internal_call <- NULL
    }

    # Wrapper for the objective function (-log_posterior)
    objective_fun <- function(params, ...) {
      res <- tryCatch(-do.call(log_posterior, c(list(params), internal_args)), error = function(e) Inf)
      if (length(res) != 1) return(Inf)
      return(res)
    }

    # Wrapper for the support function
    support_fun <- function(params, ...) {
      res <- tryCatch(do.call(log_posterior, c(list(params), internal_args)), error = function(e) -Inf)
      return(all(is.finite(res)))
    }

    if (!support_fun(x0) || !support_fun(xp0)) {
      stop("Initial points are outside the support (log-posterior is -Inf or returns an error).")
    }

    U <- objective_fun(x0); Up <- objective_fun(xp0)
    x_current <- x0; xp_current <- xp0

    x_samples <- matrix(NA, nrow = n_iter, ncol = n_dim)
    xp_samples <- matrix(NA, nrow = n_iter, ncol = n_dim)
    accepted_count <- 0

    use_progress_bar <- !is_internal_call
    if (use_progress_bar) {
      progress_bar <- utils::txtProgressBar(min = 0, max = n_iter, style = 3, width = 50, char = "=")
    }

    for (i in 1:n_iter) {
      move <- do.call(twalk_move, c(
        list(n_dim = n_dim, log_post_fun = objective_fun, support_fun = support_fun,
             x = x_current, U = U, xp = xp_current, Up = Up),
        internal_args
      ))

      if (stats::runif(1) < move$alpha) {
        x_current <- move$y
        U <- move$prop_U
        xp_current <- move$yp
        Up <- move$prop_Up
        accepted_count <- accepted_count + 1
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
    if (use_progress_bar) {
      cat(sprintf("\nAcceptance rate: %.2f%%\n", acceptance_rate * 100))
    }

    return(list(
      all_samples = rbind(x_samples, xp_samples),
      acceptance_rate = acceptance_rate,
      n_iter = n_iter,
      n_dim = n_dim
    ))
  }

  # --- PARALLEL BLOCK ---
  else {

    if (is.null(n_cores)) {
      n_cores <- max(1, parallel::detectCores() - 1)
    }
    n_cores_used <- min(n_chains, n_cores)

    cat(sprintf("--- Running t-walk in PARALLEL mode (%d chains on %d cores) ---\n", n_chains, n_cores_used))

    cl <- parallel::makeCluster(n_cores_used)
    on.exit(parallel::stopCluster(cl))

    # Export all necessary objects and functions from the global environment
    all_objects <- ls(.GlobalEnv)
    parallel::clusterExport(cl, varlist = all_objects, envir = .GlobalEnv)

    # Load required packages on each worker node
    parallel::clusterEvalQ(cl, {
      # Add any packages your log_posterior might need
      # e.g., library(mvtnorm)
    })

    # This is the function that will be executed on each worker node
    run_single_chain <- function(chain_index) {
      # Set a different seed for each chain to ensure independence
      set.seed(as.integer(Sys.time()) + chain_index)

      n_dim <- length(x0)
      # Jitter initial points slightly for each chain
      x0_i <- stats::rnorm(n_dim, mean = x0, sd = 0.1)
      xp0_i <- stats::rnorm(n_dim, mean = xp0, sd = 0.1)

      # Use 'do.call' to safely construct the function call,
      # passing the extra arguments (...) correctly.
      chain_result <- do.call(twalk, c(
        list(log_posterior = log_posterior, n_iter = n_iter, x0 = x0_i, xp0 = xp0_i,
             n_chains = 1, internal_call = TRUE),
        extra_args
      ))
      return(chain_result)
    }

    cat("Distributing work among cores...\n")
    results_list <- parallel::parLapply(cl, 1:n_chains, run_single_chain)

    cat("Chains completed. Combining results...\n")

    combined_samples <- do.call(rbind, lapply(results_list, function(res) res$all_samples))
    mean_acceptance_rate <- mean(sapply(results_list, function(res) res$acceptance_rate))
    cat(sprintf("\nMean acceptance rate across chains: %.2f%%\n", mean_acceptance_rate * 100))

    return(list(
      all_samples = combined_samples,
      acceptance_rate = mean_acceptance_rate,
      total_iterations = n_iter * n_chains,
      n_dim = length(x0),
      individual_chains = results_list
    ))
  }
}
