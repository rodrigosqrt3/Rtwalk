[![CRAN status](https://www.r-pkg.org/badges/version/Rtwalk)](https://CRAN.R-project.org/package=Rtwalk) [![R-CMD-check](https://github.com/rodrigosqrt3/Rtwalk/actions/workflows/r.yml/badge.svg)](https://github.com/rodrigosqrt3/Rtwalk/actions/workflows/r.yml) [![Codecov test coverage](https://codecov.io/github/rodrigosqrt3/Rtwalk/graph/badge.svg)](https://app.codecov.io/github/rodrigosqrt3/Rtwalk)

# Rtwalk: An MCMC Sampler Using the t-walk Algorithm

`Rtwalk` is an R implementation of the **t-walk**, a general-purpose Markov Chain Monte Carlo (MCMC) sampler for continuous distributions, ideal for Bayesian inference problems.

## Overview

This package provides an implementation of the t-walk algorithm, as originally proposed by Christen & Fox (2010). The t-walk is a robust, self-adjusting MCMC sampler, which means it does not require the tedious manual tuning of proposal parameters. It is designed to efficiently explore a wide range of target distributions, maintaining good performance even in high-dimensional or multimodal problems.

The implementation supports reproducible sequential and parallel chains,
convergence summaries, visualization helpers, constrained parameter spaces,
and numerically stable Metropolis--Hastings acceptance probabilities.

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
  log(0.5 * dnorm(x, mean = -3, sd = 0.5) + 0.5 * dnorm(x, mean = 3, sd = 0.5))
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
auxiliary t-walk trajectory is available as `result$companion_samples`.
`result$all_samples` is retained for backward compatibility, but should not be
interpreted as one time-ordered chain.

Acceptance and movement can be inspected separately:

```r
unlist(result[c("acceptance_rate", "move_rate", "no_move_rate")])
```

Here, `acceptance_rate` is the Metropolis--Hastings acceptance rate,
`move_rate` counts accepted proposals that changed the state, and
`no_move_rate` counts accepted identity proposals from the original t-walk
transition kernel.

## Parallel chains

Multiple reproducible chains can be run in parallel. Calling `set.seed()`
before `twalk()` controls the parallel random-number streams:

```r
set.seed(123)

parallel_result <- twalk(
  log_posterior = log_posterior_bimodal,
  n_iter = 25000,
  x0 = initial_point_1,
  xp0 = initial_point_2,
  n_chains = 2,
  n_cores = 2,
  show_progress = FALSE
)

summary(parallel_result, burnin_frac = 0.2)
```

## Citation

This package is an implementation of the algorithm described in the following paper. If you use `Rtwalk` in your research, please cite the original work:

> Christen, J. A., & Fox, C. (2010). A general purpose sampling algorithm for continuous distributions (the t-walk). *Bayesian Analysis*, 5(2), 263-282. [doi:10.1214/10-BA603](https://doi.org/10.1214/10-BA603)

## License

This package is licensed under the GPL-3. See the `LICENSE` file for more details.


