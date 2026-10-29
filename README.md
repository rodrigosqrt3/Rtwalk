[![CRAN status](https://www.r-pkg.org/badges/version/Rtwalk)](https://CRAN.R-project.org/package=Rtwalk) [![R-CMD-check](https://github.com/rodrigosqrt3/Rtwalk/actions/workflows/r.yml/badge.svg)](https://github.com/rodrigosqrt3/Rtwalk/actions/workflows/r.yml) [![Codecov test coverage](https://codecov.io/github/rodrigosqrt3/Rtwalk/graph/badge.svg)](https://app.codecov.io/github/rodrigosqrt3/Rtwalk)

# Rtwalk: An MCMC Sampler Using the t-walk Algorithm

`Rtwalk` is an R implementation of the **t-walk**, a general-purpose Markov Chain Monte Carlo (MCMC) sampler for continuous distributions, ideal for Bayesian inference problems.

## Overview

This package provides an implementation of the t-walk algorithm, as originally proposed by Christen & Fox (2010). The t-walk is a robust, self-adjusting MCMC sampler, which means it does not require the tedious manual tuning of proposal parameters. It is designed to efficiently explore a wide range of target distributions, maintaining good performance even in high-dimensional or multimodal problems.

Highlights:

* **Faithful to the paper.** The default settings are those of Christen & Fox
  (2010). A pure-R reference implementation follows the paper line by line,
  and the compiled C++ engine is tested to reproduce it draw by draw.
* **Fast.** The overhead per iteration is about 2 microseconds, so the cost of
  sampling is dominated by your log posterior.
* **Modern diagnostics.** Rank-normalized R-hat and bulk/tail effective sample
  sizes (Vehtari et al., 2021), integrated autocorrelation times and Monte
  Carlo standard errors, plus conversion to `coda` and `posterior` objects.
* **Practical tools.** Reproducible parallel chains with dispersed starting
  points, thinning, continuation of runs (`twalk_continue()`), single
  transitions for Gibbs hybrids (`twalk_step()`), constrained supports and
  named parameters.

## Installation

The stable version of Rtwalk is available on CRAN and can be installed with:

```r
install.packages("Rtwalk")
```

## Example Usage

Here is a simple example of how to use `Rtwalk` to sample from a bimodal distribution.

```r
library(Rtwalk)

# Define the log-posterior density of the target distribution
# In this case, a mixture of two normal distributions
log_posterior_bimodal <- function(x) {
  log_component_1 <- log(0.5) + dnorm(x, mean = -3, sd = 0.5, log = TRUE)
  log_component_2 <- log(0.5) + dnorm(x, mean = 3, sd = 0.5, log = TRUE)
  max_log <- max(log_component_1, log_component_2)

  max_log + log(
    exp(log_component_1 - max_log) + exp(log_component_2 - max_log)
  )
}

initial_point_1 <- -3
initial_point_2 <- 3

result <- twalk(
  log_posterior = log_posterior_bimodal,
  n_iter = 50000,
  x0 = initial_point_1,
  xp0 = initial_point_2
)

burnin <- floor(nrow(result$samples) * 0.2)
samples <- result$samples[
  seq.int(burnin + 1L, nrow(result$samples)),
  ,
  drop = FALSE
]

par(mfrow = c(1, 2))
hist(samples, breaks = 50, freq = FALSE, 
    main = "Posterior Distribution", xlab = "Parameter Value")
lines(density(samples), col = "blue", lwd = 2)

plot(samples, type = 'l', col = "grey30", 
    main = "Trace Plot", xlab = "Iteration")
```

The primary time-ordered trajectory is stored in `result$samples`. The
auxiliary t-walk trajectory is available as `result$companion_samples`, and
the log posterior of the stored points as `result$log_posterior`.

```r
summary(result)   # posterior summaries, ESS, MCSE, IAT and R-hat
plot(result)      # traces of the log posterior and of the parameters
```

Acceptance and movement can be inspected separately:

```r
unlist(result[c("acceptance_rate", "move_rate", "no_move_rate")])
```

Here, `acceptance_rate` is the Metropolis--Hastings acceptance rate,
`move_rate` counts accepted proposals that changed the state, and
`no_move_rate` counts accepted identity proposals from the original t-walk
transition kernel.

Proposal behavior can also be inspected separately for each t-walk kernel:

```r
result$kernel_statistics
```

This table reports the number of proposed, accepted, and state-changing moves
for the Traverse, Walk, Blow, and Hop kernels, together with their individual
acceptance and movement rates.

`print(result)` reports the overall sampler rates, while `summary(result)`
returns a table of posterior summaries and convergence diagnostics (bulk and
tail ESS, MCSE, IAT and rank-normalized R-hat), printed together with the
per-kernel performance.

## Parallel chains

Multiple reproducible chains can be run in parallel. Give one starting point
per chain (as matrix rows) so that R-hat can detect lack of convergence.
Calling `set.seed()` before `twalk()` makes the results reproducible, whatever
the number of cores:

```r
set.seed(123)

parallel_result <- twalk(
  log_posterior = log_posterior_bimodal,
  n_iter = 25000,
  x0 = matrix(c(-4, -2, 2, 4), ncol = 1),
  xp0 = matrix(c(-3, -1, 1, 3), ncol = 1),
  n_chains = 4,
  n_cores = 2,
  show_progress = FALSE
)

summary(parallel_result, burnin_frac = 0.2)
```

## Longer runs and fine control

```r
# Continue the same Markov chain for 50,000 more iterations
result <- twalk_continue(result, log_posterior_bimodal, n_iter = 50000)

# Store one iteration in ten and drop the companion trajectory
light <- twalk(log_posterior_bimodal, n_iter = 1e6, x0 = -3, xp0 = 3,
               thin = 10, control = twalk_control(store_companion = FALSE))
```

The constants of the algorithm are available through `twalk_control()`; their
defaults are those recommended by Christen & Fox (2010) and rarely need to be
changed.

## Vignettes

* `vignette("validation_study", package = "Rtwalk")`: the workflow of the
  package on targets with known answers.
* `vignette("reproducing_christen_fox", package = "Rtwalk")`: the experiments
  of the original article, reproduced with this implementation.
* `vignette("sampler_comparison", package = "Rtwalk")`: a benchmark against
  random-walk and adaptive Metropolis, an affine-invariant ensemble sampler
  and slice sampling.

Results can also be converted directly to standard `coda` objects:

```r
single_chain <- coda::as.mcmc(result, burnin_frac = 0.2)
parallel_chains <- coda::as.mcmc.list(
  parallel_result,
  burnin_frac = 0.2
)
```

The conversion preserves chain boundaries for parallel output. Set
`companion = TRUE` to convert the auxiliary trajectory instead.

## Citation

This package is an implementation of the algorithm described in the following paper. If you use `Rtwalk` in your research, please cite the original work:

> Christen, J. A., & Fox, C. (2010). A general purpose sampling algorithm for continuous distributions (the t-walk). *Bayesian Analysis*, 5(2), 263-282. [doi:10.1214/10-BA603](https://doi.org/10.1214/10-BA603)

## License

This package is licensed under GPL-3.


