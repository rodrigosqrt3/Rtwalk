# Rtwalk 2.2.0

The algorithm is unchanged: it still follows Christen & Fox (2010) exactly,
with the same default settings. Because random numbers are now drawn in a
different order, chains differ from those of version 2.1.0 for the same seed.

## Performance and memory

* The t-walk kernel is now implemented in C++ (via Rcpp). The overhead per
  iteration fell from about 65 to about 2 microseconds, so cheap targets
  sample about 30 times faster.
* Errors raised by `log_posterior` are caught without wrapping every
  evaluation in `tryCatch()`, which alone cost several times the evaluation
  of a simple target.
* The pure-R implementation is kept as a reference that follows the paper
  literally (`twalk_control(engine = "R")`); tests check that both engines
  produce the same chains draw by draw.
* New `thin` argument to store only every `thin`-th iteration.
* `twalk_control(store_companion = FALSE)` drops the companion trajectory,
  halving the memory used by the output.
* `all_samples` is no longer stored, which used as much memory as the two
  trajectories together. It is deprecated: accessing it rebuilds it, with a
  warning.

## New features

* `twalk_control()` exposes the constants of the algorithm (`n1`, `aw`, `at`
  and the move probabilities), with the defaults of Christen & Fox (2010,
  Section 2.5), for sensitivity studies.
* `twalk_continue()` extends a previous run from the final state of the
  sampler.
* `twalk_step()` performs a single t-walk transition, so the kernel can be
  combined with Gibbs updates as described in Christen & Fox (2010,
  Section 5).
* The log posterior at the stored points is returned in `log_posterior` and
  `companion_log_posterior`.
* Parameter names: if `x0` has names, `log_posterior` receives named vectors
  and the output columns are named.
* `x0` and `xp0` may be matrices with one row per chain, so that multiple
  chains start from dispersed points. A note is shown when all chains share
  the same starting points.
* Each chain now uses its own random-number stream, so multi-chain results
  are reproducible with `set.seed()` and do not depend on `n_cores`. With
  `n_cores = 1`, chains run sequentially without starting a cluster.
* Interrupting a sequential run returns the iterations completed so far.
* Errors raised by `log_posterior` at proposals are still treated as
  rejections, but a warning now reports how many occurred and the first error
  message.
* `plot()` method showing traces of the log posterior and of the parameters,
  one colour per chain.
* `posterior::as_draws()`, `as_draws_array()` and `as_draws_df()` methods.

## Diagnostics

* `calculate_diagnostics()` and `summary()` now use the rank-normalized
  split-R-hat and the bulk and tail effective sample sizes of Vehtari et al.
  (2021); the estimates agree with those of the posterior package. R-hat is
  now also available for a single chain.
* `summary()` returns a `summary.twalk` data frame, printed with the sampler
  performance, instead of printing and returning the table invisibly.
* Diagnostics use the column names of the samples, if any.

## Documentation

* New vignette reproducing the experiments of Christen & Fox (2010).
* New vignette comparing the t-walk with other general-purpose samplers.
* The validation vignette now shows its results next to the true values.

## Other changes in this version

* Refactored sequential and parallel execution to use one shared chain engine.
* Limited progress-bar refreshes to approximately 100 per chain.
* Added an explicit error for identical initial points, which otherwise create
  a degenerate t-walk state.
* Added `coda::as.mcmc()` and `coda::as.mcmc.list()` methods for single and
  parallel t-walk results, with optional burn-in and companion conversion.
* Parallel workers load the same version of the package as the main process,
  including development versions loaded with `pkgload::load_all()`.
* Added per-kernel proposal, acceptance, and actual-movement statistics for the
  Traverse, Walk, Blow, and Hop kernels.
* Corrected the Hop proposal-density helper so that only selected coordinates
  contribute to its continuous density.
* Added regression tests for Hop coordinate selection and the scalar affine
  invariance property established by Christen & Fox (2010).

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
