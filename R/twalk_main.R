.detect_available_cores <- function(detector = parallel::detectCores) {
  detected_cores <- detector()
  if (is.na(detected_cores)) {
    detected_cores <- 1L
  }

  detected_cores
}

.summarize_kernel_counts <- function(proposed, accepted, moved) {
  data.frame(
    kernel = names(proposed),
    proposed = as.integer(proposed),
    accepted = as.integer(accepted),
    moved = as.integer(moved),
    acceptance_rate = ifelse(proposed > 0L, accepted / proposed, NA_real_),
    move_rate = ifelse(proposed > 0L, moved / proposed, NA_real_),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

.is_positive_integer <- function(x) {
  is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
    x >= 1 && x == floor(x) && x <= .Machine$integer.max
}

#' Control parameters for the t-walk
#'
#' Collects the tuning constants of the t-walk. The defaults are the values
#' recommended by Christen & Fox (2010, Section 2.5), obtained by minimizing
#' integrated autocorrelation times over a wide range of targets; the authors'
#' intention is that they are left unchanged. Changing them is mainly useful
#' for sensitivity studies.
#'
#' @param n1 Expected number of coordinates moved at each iteration. Each
#'   coordinate is updated with probability `min(n, n1) / n`, where `n` is the
#'   dimension (Section 2). The paper explored `n1` in `[2, 20]`.
#' @param aw Walk-move parameter `a_w > 0` (Section 2.1). The paper explored
#'   `[0.3, 2]`.
#' @param at Traverse-move parameter `a_t > 1` (Section 2.2). The paper
#'   explored `[2, 10]`.
#' @param kernel_probs Probabilities of choosing the Traverse, Walk, Blow and
#'   Hop moves, in this order (or named). They are normalized to sum to one.
#'   The defaults correspond to the move ratios 60:60:1:1.
#' @param store_companion Logical; whether to store the trajectory of the
#'   companion point `x'`. Setting it to `FALSE` halves the memory used by
#'   the output.
#' @param engine Either `"cpp"` (default), the compiled engine, or `"R"`, a
#'   much slower pure-R reference implementation that produces identical
#'   chains and is used to validate the compiled engine.
#'
#' @return An object of class `twalk_control`.
#'
#' @references Christen, J. A. and Fox, C. (2010). A general purpose sampling
#'   algorithm for continuous distributions (the t-walk). *Bayesian Analysis*,
#'   5(2), 263--282. \doi{10.1214/10-BA603}
#'
#' @examples
#' twalk_control()
#'
#' # A sensitivity study on the Walk-move parameter
#' twalk_control(aw = 0.5)
#' @export
twalk_control <- function(n1 = 4, aw = 1.5, at = 6,
                          kernel_probs = c(traverse = 0.4918, walk = 0.4918,
                                           blow = 0.0082, hop = 0.0082),
                          store_companion = TRUE,
                          engine = c("cpp", "R")) {
  is_number <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)

  if (!is_number(n1) || n1 <= 0) {
    stop("`n1` must be one positive number.", call. = FALSE)
  }
  if (!is_number(aw) || aw <= 0) {
    stop("`aw` must be one positive number.", call. = FALSE)
  }
  if (!is_number(at) || at <= 1) {
    stop("`at` must be one number greater than 1.", call. = FALSE)
  }
  if (!is.numeric(kernel_probs) || length(kernel_probs) != 4L ||
      any(!is.finite(kernel_probs)) || any(kernel_probs < 0) ||
      sum(kernel_probs) <= 0) {
    stop(
      "`kernel_probs` must be four non-negative numbers with a positive sum.",
      call. = FALSE
    )
  }
  if (!is.null(names(kernel_probs))) {
    if (!setequal(names(kernel_probs), .kernel_names)) {
      stop(
        "Names of `kernel_probs` must be 'traverse', 'walk', 'blow' and 'hop'.",
        call. = FALSE
      )
    }
    kernel_probs <- kernel_probs[.kernel_names]
  }
  if (!is.logical(store_companion) || length(store_companion) != 1L ||
      is.na(store_companion)) {
    stop("`store_companion` must be TRUE or FALSE.", call. = FALSE)
  }

  kernel_probs <- stats::setNames(
    as.numeric(kernel_probs) / sum(kernel_probs),
    .kernel_names
  )

  structure(
    list(
      n1 = n1,
      aw = aw,
      at = at,
      kernel_probs = kernel_probs,
      store_companion = store_companion,
      engine = match.arg(engine)
    ),
    class = "twalk_control"
  )
}

#' @export
print.twalk_control <- function(x, ...) {
  cat("t-walk control parameters\n")
  cat(sprintf("  n1 = %g, aw = %g, at = %g\n", x$n1, x$aw, x$at))
  cat(
    "  kernel probabilities:",
    paste(sprintf("%s = %.4f", names(x$kernel_probs), x$kernel_probs),
          collapse = ", "),
    "\n"
  )
  cat("  store companion:", x$store_companion, "\n")
  cat("  engine:", x$engine, "\n")
  invisible(x)
}

.validate_control <- function(control) {
  if (!inherits(control, "twalk_control")) {
    stop("`control` must be created with `twalk_control()`.", call. = FALSE)
  }
  control
}

# Build f(params) = log_posterior(params, <extra arguments>) once, so that
# each evaluation is a plain closure call.
.make_log_posterior <- function(log_posterior, extra_args) {
  if (length(extra_args) == 0L) {
    return(log_posterior)
  }

  bind_arguments <- function(...) {
    force(log_posterior)
    function(params) log_posterior(params, ...)
  }
  do.call(bind_arguments, extra_args)
}

.initial_log_posterior <- function(eval_fun, x) {
  value <- tryCatch(eval_fun(x), error = function(e) e)
  if (inherits(value, "error")) {
    return(list(value = -Inf, error = conditionMessage(value)))
  }
  list(value = .check_log_posterior_value(value), error = NULL)
}

# Run one chain from (x0, xp0). `lp0` holds the log posterior at the initial
# points when already known (continuation and parallel runs).
.run_chain <- function(eval_fun, x0, xp0, n_iter, thin, control,
                       show_progress = FALSE, lp0 = NULL,
                       progress_label = NULL) {
  n_dim <- length(x0)
  param_names <- names(x0)
  names(xp0) <- param_names

  if (is.null(lp0)) {
    init_x <- .initial_log_posterior(eval_fun, x0)
    init_xp <- .initial_log_posterior(eval_fun, xp0)
    lp0 <- c(init_x$value, init_xp$value)

    if (!all(is.finite(lp0))) {
      init_error <- c(init_x$error, init_xp$error)
      stop(
        "Initial points are outside the support (log-posterior is -Inf or returns an error).",
        if (length(init_error) > 0L) paste0("\nFirst error: ", init_error[[1L]]),
        call. = FALSE
      )
    }
  }

  engine_new <- if (identical(control$engine, "R")) .engine_r_new else .engine_cpp_new
  engine <- engine_new(
    x0 = x0, xp0 = xp0, lp_x0 = lp0[[1L]], lp_xp0 = lp0[[2L]],
    n_iter = n_iter, thin = thin,
    p_phi = min(n_dim, control$n1) / n_dim,
    aw = control$aw, at = control$at,
    cum_probs = cumsum(control$kernel_probs),
    store_companion = control$store_companion,
    names = param_names
  )

  if (show_progress) {
    if (!is.null(progress_label)) message(progress_label)
    progress_bar <- utils::txtProgressBar(
      min = 0, max = n_iter, style = 3, width = 50, char = "="
    )
    on.exit(close(progress_bar), add = TRUE)
    chunk <- max(1L, as.integer(ceiling(n_iter / 100)))
  } else {
    chunk <- n_iter
  }

  n_errors <- 0L
  first_error <- NULL
  interrupted <- FALSE

  while (engine$iter() < n_iter) {
    outcome <- tryCatch(
      {
        engine$run(eval_fun, chunk)
        NULL
      },
      error = function(e) e,
      interrupt = function(e) e
    )

    if (inherits(outcome, "interrupt")) {
      interrupted <- TRUE
      break
    }

    if (inherits(outcome, "error")) {
      if (!engine$pending()) {
        stop(conditionMessage(outcome), call. = FALSE)
      }
      n_errors <- n_errors + 1L
      if (is.null(first_error)) first_error <- conditionMessage(outcome)
      engine$reject_pending()
    }

    if (show_progress) {
      utils::setTxtProgressBar(progress_bar, engine$iter())
    }
  }

  res <- engine$result()
  iterations <- res$iter

  if (!is.null(param_names)) {
    colnames(res$samples) <- param_names
    if (!is.null(res$companion_samples)) colnames(res$companion_samples) <- param_names
    names(res$x) <- param_names
    names(res$xp) <- param_names
  }

  proposed <- stats::setNames(res$proposed, .kernel_names)
  accepted <- stats::setNames(res$accepted, .kernel_names)
  moved <- stats::setNames(res$moved, .kernel_names)
  denominator <- max(iterations, 1L)

  structure(
    list(
      samples = res$samples,
      companion_samples = res$companion_samples,
      log_posterior = res$log_posterior,
      companion_log_posterior = res$companion_log_posterior,
      acceptance_rate = sum(accepted) / denominator,
      move_rate = sum(moved) / denominator,
      no_move_rate = (sum(accepted) - sum(moved)) / denominator,
      kernel_statistics = .summarize_kernel_counts(proposed, accepted, moved),
      n_iter = iterations,
      thin = thin,
      n_chains = 1L,
      total_iterations = iterations,
      n_dim = n_dim,
      n_evaluations = res$n_evals,
      n_errors = n_errors,
      first_error = first_error,
      interrupted = interrupted,
      final_state = list(
        x = res$x, xp = res$xp,
        log_x = res$lp_x, log_xp = res$lp_xp
      ),
      control = control
    ),
    class = "twalk"
  )
}

# Turn x0 / xp0 (a vector shared by all chains, or a matrix with one row per
# chain) into a list of per-chain starting points.
.chain_starts <- function(x, n_chains, arg) {
  if (is.matrix(x)) {
    if (!is.numeric(x) || nrow(x) != n_chains || ncol(x) < 1L ||
        any(!is.finite(x))) {
      stop(
        sprintf(
          "When `%s` is a matrix it must be numeric, finite and have one row per chain (%d).",
          arg, n_chains
        ),
        call. = FALSE
      )
    }
    return(lapply(seq_len(n_chains), function(i) {
      stats::setNames(x[i, ], colnames(x))
    }))
  }

  if (!is.numeric(x) || !is.null(dim(x)) || length(x) == 0L ||
      any(!is.finite(x))) {
    stop(
      "`x0` and `xp0` must be non-empty numeric vectors with finite values, or matrices with one row per chain.",
      call. = FALSE
    )
  }

  rep(list(x), n_chains)
}

# One independent L'Ecuyer-CMRG stream per chain, derived from the user's
# current random-number state. Chain i always receives the same stream,
# whatever the number of cores, so results do not depend on `n_cores`.
.chain_seeds <- function(n_chains) {
  # Drawing the base seed also guarantees that `.Random.seed` exists.
  base_seed <- sample.int(.Machine$integer.max, 1L)

  old_seed <- get(".Random.seed", envir = globalenv())
  old_kind <- RNGkind()
  on.exit({
    RNGkind(old_kind[[1L]], old_kind[[2L]], old_kind[[3L]])
    assign(".Random.seed", old_seed, envir = globalenv())
  })

  RNGkind("L'Ecuyer-CMRG")
  set.seed(base_seed)
  seed <- get(".Random.seed", envir = globalenv())

  seeds <- vector("list", n_chains)
  for (i in seq_len(n_chains)) {
    seeds[[i]] <- seed
    seed <- parallel::nextRNGStream(seed)
  }
  seeds
}

.with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv())
  old_kind <- RNGkind()
  on.exit({
    RNGkind(old_kind[[1L]], old_kind[[2L]], old_kind[[3L]])
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })

  RNGkind("L'Ecuyer-CMRG")
  assign(".Random.seed", seed, envir = globalenv())
  expr
}

# Run one chain from a chain specification. Used both in the main process and
# on parallel workers.
.run_chain_spec <- function(spec, eval_fun, n_iter, thin, control) {
  .with_seed(
    spec$seed,
    .run_chain(
      eval_fun = eval_fun,
      x0 = spec$x0,
      xp0 = spec$xp0,
      n_iter = n_iter,
      thin = thin,
      control = control,
      show_progress = FALSE,
      lp0 = spec$lp0
    )
  )
}

# Functions created in the global or a local environment may refer to data
# and helper functions stored alongside them. PSOCK workers start with a
# clean workspace, so recursively collect those dependencies.
.collect_dependencies <- function(fun, collected = list()) {
  fun_env <- environment(fun)
  if (is.null(fun_env) || isNamespace(fun_env) ||
      identical(fun_env, baseenv())) {
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
      collected <- .collect_dependencies(object, collected)
    }
  }

  collected
}

# Path of the package source when it was loaded with pkgload::load_all();
# NULL for an installed package. Workers must load the same version.
# These two functions only matter during development and run on workers,
# outside the reach of coverage tools.
# nocov start
.dev_package_path <- function() {
  path <- getNamespaceInfo("Rtwalk", "path")
  is_source <- file.exists(file.path(path, "DESCRIPTION")) &&
    !file.exists(file.path(path, "Meta", "package.rds"))

  if (is_source && isNamespaceLoaded("pkgload")) path else NULL
}

.setup_worker <- function(dev_path) {
  if (is.null(dev_path)) {
    loadNamespace("Rtwalk")
  } else {
    pkgload::load_all(dev_path, quiet = TRUE, export_all = FALSE,
                      helpers = FALSE, attach_testthat = FALSE)
  }
  invisible(TRUE)
}
# nocov end

.run_chains_parallel <- function(specs, log_posterior, extra_args,
                                 n_iter, thin, control, n_cores) {
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  setup <- .setup_worker
  environment(setup) <- globalenv()
  parallel::clusterCall(cl, setup, .dev_package_path())

  dependencies <- .collect_dependencies(log_posterior)
  if (length(dependencies) > 0L) {
    parallel::clusterExport(
      cl,
      varlist = names(dependencies),
      envir = list2env(dependencies, parent = emptyenv())
    )
  }

  worker <- function(spec, log_posterior, extra_args, n_iter, thin, control) {
    ns <- asNamespace("Rtwalk")
    eval_fun <- ns$.make_log_posterior(log_posterior, extra_args)
    ns$.run_chain_spec(spec, eval_fun, n_iter, thin, control)
  }
  environment(worker) <- globalenv()

  parallel::parLapply(
    cl, specs, worker,
    log_posterior = log_posterior, extra_args = extra_args,
    n_iter = n_iter, thin = thin, control = control
  )
}

.combine_chains <- function(results) {
  stack <- function(field) {
    parts <- lapply(results, `[[`, field)
    if (any(vapply(parts, is.null, logical(1)))) return(NULL)
    if (is.matrix(parts[[1L]])) do.call(rbind, parts) else unlist(parts, use.names = FALSE)
  }
  sum_counts <- function(field) {
    Reduce(`+`, lapply(results, function(res) res$kernel_statistics[[field]]))
  }

  iterations <- vapply(results, `[[`, numeric(1), "n_iter")
  errors <- unlist(lapply(results, `[[`, "first_error"))

  structure(
    list(
      samples = stack("samples"),
      companion_samples = stack("companion_samples"),
      log_posterior = stack("log_posterior"),
      companion_log_posterior = stack("companion_log_posterior"),
      acceptance_rate = mean(vapply(results, `[[`, numeric(1), "acceptance_rate")),
      move_rate = mean(vapply(results, `[[`, numeric(1), "move_rate")),
      no_move_rate = mean(vapply(results, `[[`, numeric(1), "no_move_rate")),
      kernel_statistics = .summarize_kernel_counts(
        stats::setNames(sum_counts("proposed"), .kernel_names),
        stats::setNames(sum_counts("accepted"), .kernel_names),
        stats::setNames(sum_counts("moved"), .kernel_names)
      ),
      n_iter = results[[1L]]$n_iter,
      thin = results[[1L]]$thin,
      n_chains = length(results),
      total_iterations = sum(iterations),
      n_dim = results[[1L]]$n_dim,
      n_evaluations = sum(vapply(results, `[[`, numeric(1), "n_evaluations")),
      n_errors = sum(vapply(results, `[[`, numeric(1), "n_errors")),
      first_error = if (length(errors) > 0L) errors[[1L]],
      interrupted = any(vapply(results, `[[`, logical(1), "interrupted")),
      control = results[[1L]]$control,
      individual_chains = results
    ),
    class = "twalk"
  )
}

.warn_sampling_problems <- function(result) {
  if (isTRUE(result$n_errors > 0)) {
    warning(
      sprintf(
        paste0(
          "`log_posterior` raised an error in %d of %d evaluations; those ",
          "proposals were treated as outside the support and rejected. ",
          "Return -Inf outside the support to silence this warning. ",
          "First error: %s"
        ),
        result$n_errors, result$n_evaluations + result$n_errors,
        result$first_error
      ),
      call. = FALSE
    )
  }
  if (isTRUE(result$interrupted)) {
    warning(
      sprintf(
        "Sampling was interrupted; returning the %d completed iterations.",
        result$total_iterations
      ),
      call. = FALSE
    )
  }
}

# Shared driver for twalk() and twalk_continue().
.twalk_run <- function(log_posterior, extra_args, starts, n_iter, thin,
                       control, n_chains, n_cores, show_progress,
                       chain_labels = TRUE) {
  eval_fun <- .make_log_posterior(log_posterior, extra_args)

  if (n_chains == 1L) {
    spec <- starts[[1L]]
    if (show_progress) {
      message("--- Running t-walk in sequential mode (1 chain) ---")
    }
    result <- .run_chain(
      eval_fun, spec$x0, spec$xp0, n_iter, thin, control,
      show_progress = show_progress, lp0 = spec$lp0
    )
    if (show_progress) {
      message(sprintf("\nAcceptance rate: %.2f%%", result$acceptance_rate * 100))
      message(sprintf("Actual move rate: %.2f%%", result$move_rate * 100))
    }
    .warn_sampling_problems(result)
    return(result)
  }

  # Evaluate starting points in the main process so that problems are
  # reported before any worker starts.
  for (i in seq_len(n_chains)) {
    if (is.null(starts[[i]]$lp0)) {
      init_x <- .initial_log_posterior(eval_fun, starts[[i]]$x0)
      init_xp <- .initial_log_posterior(eval_fun, starts[[i]]$xp0)
      lp0 <- c(init_x$value, init_xp$value)
      if (!all(is.finite(lp0))) {
        init_error <- c(init_x$error, init_xp$error)
        stop(
          sprintf("Initial points of chain %d are outside the support (log-posterior is -Inf or returns an error).", i),
          if (length(init_error) > 0L) paste0("\nFirst error: ", init_error[[1L]]),
          call. = FALSE
        )
      }
      starts[[i]]$lp0 <- lp0
    }
  }

  seeds <- .chain_seeds(n_chains)
  specs <- lapply(seq_len(n_chains), function(i) {
    c(starts[[i]], list(seed = seeds[[i]]))
  })

  if (is.null(n_cores)) {
    n_cores <- max(1L, .detect_available_cores() - 1L)
  }
  n_cores_used <- min(n_chains, n_cores)

  if (n_cores_used == 1L) {
    if (show_progress) {
      message(sprintf("--- Running t-walk sequentially (%d chains on 1 core) ---", n_chains))
    }
    results <- lapply(seq_len(n_chains), function(i) {
      spec <- specs[[i]]
      .with_seed(
        spec$seed,
        .run_chain(
          eval_fun, spec$x0, spec$xp0, n_iter, thin, control,
          show_progress = show_progress, lp0 = spec$lp0,
          progress_label = if (show_progress) sprintf("Chain %d of %d", i, n_chains)
        )
      )
    })
  } else {
    if (show_progress) {
      message(sprintf(
        "--- Running t-walk in PARALLEL mode (%d chains on %d cores) ---",
        n_chains, n_cores_used
      ))
      message("Distributing work among cores...")
    }
    results <- .run_chains_parallel(
      specs, log_posterior, extra_args, n_iter, thin, control, n_cores_used
    )
    if (show_progress) {
      message("Chains completed. Combining results...")
    }
  }

  combined <- .combine_chains(results)
  if (show_progress) {
    message(sprintf("\nMean acceptance rate across chains: %.2f%%", combined$acceptance_rate * 100))
    message(sprintf("Mean actual move rate across chains: %.2f%%", combined$move_rate * 100))
  }
  .warn_sampling_problems(combined)
  combined
}

#' Run the t-walk MCMC Algorithm
#'
#' This function implements the t-walk algorithm by Christen & Fox (2010),
#' a general-purpose MCMC sampler that does not require manual tuning.
#' The function can run multiple independent MCMC chains in parallel
#' to accelerate execution and facilitate convergence diagnostics.
#'
#' The t-walk maintains two points, `x` and `x'`, and is a Metropolis--Hastings
#' sampler on the product space with target \eqn{\pi(x)\pi(x')}. At each
#' iteration one of four moves (Traverse, Walk, Blow, Hop) is chosen and only
#' one of the two points is updated. The sampler is invariant to changes of
#' scale and location of the target, so no tuning is required. The kernel is
#' implemented in C++; [twalk_control()] gives access to the constants of the
#' algorithm and to a pure-R reference implementation.
#'
#' @section Errors in `log_posterior`:
#' Return `-Inf` outside the support of the target. If `log_posterior` raises
#' an error at a proposal, the proposal is rejected (as if outside the
#' support) and a single warning summarizing the errors is issued at the end.
#'
#' @section Interrupting:
#' In sequential mode, interrupting the computation (for example with Esc or
#' Ctrl+C) returns the iterations completed so far, with a warning.
#'
#' @param log_posterior A function that takes a parameter vector as its
#'   first argument and returns one numeric log posterior density. It may
#'   return `-Inf` outside the support, but must not return vectors, `NA`,
#'   `NaN`, or `+Inf`. If `x0` has names, the parameter vector passed to
#'   `log_posterior` carries the same names.
#'   Additional arguments can be passed to this function via `...`.
#' @param n_iter The number of iterations to run for each chain.
#' @param x0 A numeric vector with the initial values for the first point
#'   (`x`), or a matrix with one row per chain.
#' @param xp0 A numeric vector with the initial values for the second point
#'   (`x'`), or a matrix with one row per chain. It must differ from `x0`;
#'   identical initial points produce a degenerate t-walk state.
#' @param n_chains The number of independent MCMC chains to run.
#'   Defaults to `1`, which runs a single chain sequentially. If greater
#'   than 1, parallel mode is activated. Convergence diagnostics such as
#'   R-hat are most informative when the chains start from dispersed points,
#'   supplied as matrices in `x0` and `xp0`.
#' @param n_cores The number of CPU cores to use in parallel mode.
#'   If `NULL` (default), it will attempt to use all available cores minus one.
#'   Each chain receives its own random-number stream derived from R's current
#'   random state, so calling `set.seed()` before `twalk()` makes results
#'   reproducible, and the results do not depend on `n_cores`.
#' @param show_progress Logical; whether to display progress bars and status
#'   messages. Defaults to `TRUE`.
#' @param thin Positive integer; only every `thin`-th iteration is stored.
#'   All iterations are still run.
#' @param control A list of algorithm constants created by [twalk_control()].
#'   The defaults are those recommended by Christen & Fox (2010).
#' @param ... Additional arguments to be passed to the `log_posterior` function.
#'
#' @return An object of class `twalk`, a list containing:
#' \item{samples}{The primary t-walk trajectory, with `n_iter / thin` rows per
#'       chain. Columns are named after `x0`, if it has names.}
#' \item{companion_samples}{The auxiliary trajectory maintained by the t-walk,
#'       or `NULL` if `store_companion = FALSE` in [twalk_control()].}
#' \item{log_posterior, companion_log_posterior}{The log posterior evaluated at
#'       the stored points of each trajectory, useful to monitor burn-in.}
#' \item{acceptance_rate}{The average Metropolis--Hastings acceptance rate
#'       across all chains. Accepted identity proposals are included, as
#'       required by the t-walk transition kernel.}
#' \item{move_rate}{The proportion of iterations in which an accepted
#'       proposal actually changed at least one of the two t-walk points.}
#' \item{no_move_rate}{The proportion of iterations containing an accepted
#'       identity proposal. It equals `acceptance_rate - move_rate`.}
#' \item{kernel_statistics}{A data frame reporting the number of proposals,
#'       acceptances, actual movements, acceptance rate, and movement rate for
#'       each of the Traverse, Walk, Blow, and Hop kernels.}
#' \item{n_iter}{The number of iterations run per chain.}
#' \item{thin}{The thinning interval.}
#' \item{n_chains}{The number of independent chains.}
#' \item{total_iterations}{The total number of iterations run over all chains
#'       (`n_iter * n_chains`).}
#' \item{n_dim}{The dimension of the parameter space.}
#' \item{n_evaluations}{The number of successful log-posterior evaluations.}
#' \item{n_errors, first_error}{The number of evaluations that raised an error
#'       and the first error message.}
#' \item{final_state}{For single-chain results, the final points and their
#'       log posterior values, used by [twalk_continue()].}
#' \item{control}{The [twalk_control()] settings used.}
#' \item{individual_chains}{If `n_chains > 1`, a list containing the raw
#'       results from each separate chain, useful for diagnostics like R-hat.}
#'
#' The field `all_samples` of earlier versions is deprecated: it is no longer
#' stored, and accessing it computes it on the fly with a warning.
#'
#' @seealso [twalk_continue()] to extend a run, [twalk_step()] to use the
#'   t-walk kernel inside other samplers, [summary.twalk()], [plot.twalk()].
#'
#' @references Christen, J. A. and Fox, C. (2010). A general purpose sampling
#'   algorithm for continuous distributions (the t-walk). *Bayesian Analysis*,
#'   5(2), 263--282. \doi{10.1214/10-BA603}
#'
#' @export
#' @useDynLib Rtwalk, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom parallel detectCores makeCluster clusterExport clusterCall parLapply stopCluster nextRNGStream
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
#'   # Set a seed for reproducibility
#'   set.seed(123)
#'   result_seq <- twalk(log_posterior = log_post, n_iter = 5000,
#'                       x0 = c(-1, 1), xp0 = c(1, -1))
#'
#'   plot(result_seq$samples, pch = '.', main = "t-walk Samples (Sequential)")
#' }
#'
#' # Example 2: named parameters and extra arguments
#' log_post_normal <- function(theta, y) {
#'   sum(dnorm(y, theta[["mu"]], exp(theta[["log_sigma"]]), log = TRUE))
#' }
#' y <- c(4.1, 5.3, 3.8, 4.9, 5.6)
#' fit <- twalk(log_post_normal, n_iter = 5000,
#'              x0 = c(mu = 4, log_sigma = 0), xp0 = c(mu = 5, log_sigma = 1),
#'              y = y, show_progress = FALSE)
#' summary(fit)
#'
#' \donttest{
#' # Example 3: two chains in parallel from dispersed starting points
#' if (requireNamespace("mvtnorm", quietly = TRUE)) {
#'   set.seed(123)
#'   result_par <- twalk(log_posterior = log_post, n_iter = 2500,
#'                       x0 = rbind(c(-3, 3), c(3, -3)),
#'                       xp0 = rbind(c(-2, 2), c(2, -2)),
#'                       n_chains = 2, n_cores = 2)
#'   summary(result_par)
#' }
#' }
twalk <- function(log_posterior, n_iter, x0, xp0,
                  n_chains = 1, n_cores = NULL, show_progress = TRUE,
                  thin = 1, control = twalk_control(), ...) {

  if (!is.function(log_posterior)) {
    stop("`log_posterior` must be a function.", call. = FALSE)
  }

  if (!.is_positive_integer(n_iter)) {
    stop("`n_iter` must be a positive integer.", call. = FALSE)
  }

  if (!.is_positive_integer(n_chains)) {
    stop("`n_chains` must be a positive integer.", call. = FALSE)
  }

  if (!is.null(n_cores) && !.is_positive_integer(n_cores)) {
    stop("`n_cores` must be NULL or a positive integer.", call. = FALSE)
  }

  if (!is.logical(show_progress) || length(show_progress) != 1L || is.na(show_progress)) {
    stop("`show_progress` must be TRUE or FALSE.", call. = FALSE)
  }

  if (!.is_positive_integer(thin) || thin > n_iter) {
    stop("`thin` must be a positive integer not larger than `n_iter`.", call. = FALSE)
  }

  control <- .validate_control(control)
  n_iter <- as.integer(n_iter)
  n_chains <- as.integer(n_chains)
  thin <- as.integer(thin)
  if (!is.null(n_cores)) {
    n_cores <- as.integer(n_cores)
  }

  if (n_chains == 1L && (is.matrix(x0) || is.matrix(xp0))) {
    x0 <- if (is.matrix(x0) && nrow(x0) == 1L) stats::setNames(x0[1L, ], colnames(x0)) else x0
    xp0 <- if (is.matrix(xp0) && nrow(xp0) == 1L) stats::setNames(xp0[1L, ], colnames(xp0)) else xp0
  }
  x_starts <- .chain_starts(x0, n_chains, "x0")
  xp_starts <- .chain_starts(xp0, n_chains, "xp0")

  starts <- lapply(seq_len(n_chains), function(i) {
    x <- x_starts[[i]]
    xp <- xp_starts[[i]]
    if (length(x) != length(xp)) {
      stop("`x0` and `xp0` must have the same length.", call. = FALSE)
    }
    if (all(x == xp)) {
      stop("`x0` and `xp0` must be different initial points.", call. = FALSE)
    }
    if (is.null(names(x)) && !is.null(names(xp))) names(x) <- names(xp)
    list(x0 = x, xp0 = xp, lp0 = NULL)
  })

  if (show_progress && n_chains > 1L && !is.matrix(x0) && !is.matrix(xp0)) {
    message(
      "Note: all chains start from the same points. Supply matrices in `x0` ",
      "and `xp0` (one row per chain) with dispersed starting points to make ",
      "R-hat more informative."
    )
  }

  .twalk_run(
    log_posterior = log_posterior,
    extra_args = list(...),
    starts = starts,
    n_iter = n_iter,
    thin = thin,
    control = control,
    n_chains = n_chains,
    n_cores = n_cores,
    show_progress = show_progress
  )
}

.final_states <- function(object) {
  chains <- if (!is.null(object$individual_chains)) {
    object$individual_chains
  } else {
    list(object)
  }

  lapply(chains, function(chain) {
    state <- chain$final_state
    if (is.null(state)) {
      stop(
        "`object` does not contain the final state of the sampler; it was created by an older version of Rtwalk.",
        call. = FALSE
      )
    }
    list(
      x0 = state$x,
      xp0 = state$xp,
      lp0 = c(state$log_x, state$log_xp)
    )
  })
}

.append_chain <- function(old, new) {
  bind <- function(a, b) {
    if (is.null(a) || is.null(b)) return(NULL)
    if (is.matrix(a)) rbind(a, b) else c(a, b)
  }

  iterations <- old$n_iter + new$n_iter
  counts <- function(field) {
    old$kernel_statistics[[field]] + new$kernel_statistics[[field]]
  }
  accepted <- stats::setNames(counts("accepted"), .kernel_names)
  moved <- stats::setNames(counts("moved"), .kernel_names)

  new$samples <- bind(old$samples, new$samples)
  new$companion_samples <- bind(old$companion_samples, new$companion_samples)
  new$log_posterior <- bind(old$log_posterior, new$log_posterior)
  new$companion_log_posterior <- bind(old$companion_log_posterior, new$companion_log_posterior)
  new$kernel_statistics <- .summarize_kernel_counts(
    stats::setNames(counts("proposed"), .kernel_names), accepted, moved
  )
  new$acceptance_rate <- sum(accepted) / iterations
  new$move_rate <- sum(moved) / iterations
  new$no_move_rate <- (sum(accepted) - sum(moved)) / iterations
  new$n_iter <- iterations
  new$total_iterations <- iterations
  new$n_evaluations <- old$n_evaluations + new$n_evaluations
  new$n_errors <- old$n_errors + new$n_errors
  new
}

#' Continue a t-walk run
#'
#' Runs further iterations of a previous [twalk()] result, starting each chain
#' from the final state of the sampler. The Markov chain is continued exactly:
#' no burn-in or re-initialization takes place.
#'
#' @param object An object of class `twalk`.
#' @param log_posterior The log-posterior function used to create `object`.
#' @param n_iter Number of additional iterations per chain.
#' @param n_cores Number of cores for multi-chain results (see [twalk()]).
#' @param show_progress Logical; whether to display progress information.
#' @param append Logical; if `TRUE` (default), the new iterations are appended
#'   to those in `object`. If `FALSE`, only the new iterations are returned.
#' @param ... Additional arguments passed to `log_posterior`.
#'
#' @return An object of class `twalk`.
#'
#' @examples
#' log_post <- function(x) -0.5 * sum(x^2)
#' set.seed(1)
#' fit <- twalk(log_post, n_iter = 2000, x0 = c(-1, 1), xp0 = c(1, -1),
#'              show_progress = FALSE)
#' fit <- twalk_continue(fit, log_post, n_iter = 3000, show_progress = FALSE)
#' fit$n_iter
#' @export
twalk_continue <- function(object, log_posterior, n_iter, n_cores = NULL,
                           show_progress = TRUE, append = TRUE, ...) {
  if (!inherits(object, "twalk")) {
    stop("`object` must be a `twalk` result.", call. = FALSE)
  }
  if (!is.function(log_posterior)) {
    stop("`log_posterior` must be a function.", call. = FALSE)
  }
  if (!.is_positive_integer(n_iter)) {
    stop("`n_iter` must be a positive integer.", call. = FALSE)
  }
  if (!is.null(n_cores) && !.is_positive_integer(n_cores)) {
    stop("`n_cores` must be NULL or a positive integer.", call. = FALSE)
  }
  if (!is.logical(append) || length(append) != 1L || is.na(append)) {
    stop("`append` must be TRUE or FALSE.", call. = FALSE)
  }

  thin <- if (is.null(object$thin)) 1L else as.integer(object$thin)
  control <- if (is.null(object$control)) twalk_control() else object$control
  if (n_iter < thin) {
    stop("`n_iter` must be at least the thinning interval of `object`.", call. = FALSE)
  }

  starts <- .final_states(object)
  n_chains <- length(starts)

  new <- .twalk_run(
    log_posterior = log_posterior,
    extra_args = list(...),
    starts = starts,
    n_iter = as.integer(n_iter),
    thin = thin,
    control = control,
    n_chains = n_chains,
    n_cores = if (is.null(n_cores)) NULL else as.integer(n_cores),
    show_progress = show_progress
  )

  if (!append) {
    return(new)
  }

  if (n_chains == 1L) {
    return(.append_chain(object, new))
  }

  chains <- Map(.append_chain, object$individual_chains, new$individual_chains)
  .combine_chains(chains)
}

#' Perform a single t-walk transition
#'
#' Applies one transition of the t-walk kernel to the pair `(x, xp)`. This is
#' a building block for hybrid samplers: Christen & Fox (2010, Section 5)
#' combine the t-walk kernel with Gibbs updates for blocks of parameters whose
#' full conditionals are available.
#'
#' Since the t-walk targets \eqn{\pi(x)\pi(x')}, the log densities of both
#' points are needed. Pass them through `log_x` and `log_xp` to avoid
#' re-evaluating them; when the target changes between calls (as in a Gibbs
#' sampler, where the conditional depends on the other blocks), leave them as
#' `NULL` so they are recomputed.
#'
#' @param log_posterior Log-density function, as in [twalk()].
#' @param x,xp The current two points.
#' @param log_x,log_xp Optional log densities at `x` and `xp`.
#' @param control A [twalk_control()] object.
#' @param ... Additional arguments passed to `log_posterior`.
#'
#' @return A list with the new points `x` and `xp`, their log densities
#'   `log_x` and `log_xp`, the `kernel` used, and the logical flags
#'   `accepted` and `moved`. If `log_posterior` raised an error at the
#'   proposal, the proposal is rejected and the message is returned in
#'   `error`.
#'
#' @examples
#' # Gibbs sampler for (mu, tau) in a normal model, where the t-walk updates
#' # mu given tau and tau is drawn from its full conditional.
#' set.seed(1)
#' y <- rnorm(30, 2, 1.5)
#' log_cond_mu <- function(mu, tau) sum(dnorm(y, mu, 1 / sqrt(tau), log = TRUE))
#'
#' mu <- 0; mu_p <- 1; tau <- 1
#' draws <- matrix(NA, 2000, 2, dimnames = list(NULL, c("mu", "tau")))
#' for (i in 1:2000) {
#'   step <- twalk_step(log_cond_mu, mu, mu_p, tau = tau)
#'   mu <- step$x; mu_p <- step$xp
#'   tau <- rgamma(1, 1 + length(y) / 2, 1 + sum((y - mu)^2) / 2)
#'   draws[i, ] <- c(mu, tau)
#' }
#' colMeans(draws[-(1:500), ])
#' @export
twalk_step <- function(log_posterior, x, xp, log_x = NULL, log_xp = NULL,
                       control = twalk_control(), ...) {
  if (!is.function(log_posterior)) {
    stop("`log_posterior` must be a function.", call. = FALSE)
  }
  control <- .validate_control(control)
  if (!is.numeric(x) || !is.numeric(xp) || length(x) != length(xp) ||
      length(x) == 0L || any(!is.finite(x)) || any(!is.finite(xp))) {
    stop("`x` and `xp` must be finite numeric vectors of the same length.", call. = FALSE)
  }

  eval_fun <- .make_log_posterior(log_posterior, list(...))
  if (is.null(log_x)) log_x <- .check_log_posterior_value(eval_fun(x))
  if (is.null(log_xp)) log_xp <- .check_log_posterior_value(eval_fun(xp))
  if (!is.finite(log_x) || !is.finite(log_xp)) {
    stop("`x` and `xp` must be inside the support of the target.", call. = FALSE)
  }

  control$store_companion <- TRUE
  res <- .run_chain(
    eval_fun, x, xp, n_iter = 1L, thin = 1L, control = control,
    lp0 = c(log_x, log_xp)
  )
  stats <- res$kernel_statistics

  list(
    x = res$final_state$x,
    xp = res$final_state$xp,
    log_x = res$final_state$log_x,
    log_xp = res$final_state$log_xp,
    kernel = stats$kernel[stats$proposed == 1L],
    accepted = sum(stats$accepted) == 1L,
    moved = sum(stats$moved) == 1L,
    error = res$first_error
  )
}
