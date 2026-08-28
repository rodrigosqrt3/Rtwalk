# Rtwalk 2.1.0

## Algorithm correctness

* Corrected the Blow-kernel proposal density to include only coordinates
  selected for movement, matching Christen & Fox (2010).
* Bounded all Metropolis--Hastings acceptance probabilities to `[0, 1]` and
  moved their calculation to the log scale for numerical stability.
* Preserved accepted identity proposals from the original t-walk transition
  kernel while reporting them separately from actual state movements.
* Avoided repeated evaluations of the log posterior, reducing sequential
  execution to one evaluation per proposal plus the two initial evaluations.

## Sampling and parallel execution

* Added reproducible independent random-number streams for parallel chains;
  `set.seed()` now controls both sequential and parallel runs.
* Parallel chains now preserve the supplied initial points and support
  constrained parameter spaces correctly.
* Parallel execution now exports direct and recursively referenced
  log-posterior dependencies to worker processes.
* Added consistent `n_chains` and `total_iterations` metadata in sequential and
  parallel results.

## Output and diagnostics

* Added separate `samples` and `companion_samples` trajectories. The legacy
  `all_samples` field remains available for compatibility but is not treated as
  one time-ordered chain.
* Added `move_rate` and `no_move_rate` alongside the Metropolis--Hastings
  `acceptance_rate`.
* Applied burn-in independently to each parallel chain and computed effective
  sample sizes without joining chain boundaries.
* Unified sample and burn-in validation across diagnostics and visualization.

## Interface and validation

* Added the public `show_progress` argument. Additional arguments, including an
  argument named `internal_call`, are now passed exclusively to
  `log_posterior`.
* Added clear validation for iteration counts, chains, cores, initial points,
  log-posterior return values, burn-in fractions, and parameter names.
* Expanded regression tests for kernel mathematics, parallel reproducibility,
  constrained supports, dependency export, diagnostics, and edge cases.

# Rtwalk 2.0.2

* Added a comprehensive unit test suite using `testthat`, achieving 100% test coverage.
* Added `testthat` and `covr` to Suggested dependencies.

# Rtwalk 2.0.1

* Added S3 methods for `print()` and `summary()` for t-walk output objects.
* Improved diagnostics reporting and user-facing summaries.
* Minor internal refactoring with no changes to core algorithm behavior.

# Rtwalk 2.0.0

* Reintroduced the t-walk MCMC algorithm to CRAN.
* Modernized the codebase to comply with current CRAN policies.
* Updated documentation and examples.
* Ensured compatibility with recent versions of R.
