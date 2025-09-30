# Rtwalk: An MCMC Sampler Using the t-walk Algorithm

`Rtwalk` is an R implementation of the **t-walk**, a general-purpose Markov Chain Monte Carlo (MCMC) sampler for continuous distributions, ideal for Bayesian inference problems.

## Overview

This package provides a implementation of the t-walk algorithm, as originally proposed by Christen & Fox (2010). The t-walk is a robust, self-adjusting MCMC sampler, which means it does not require the tedious manual tuning of proposal parameters. It is designed to efficiently explore a wide range of target distributions, maintaining good performance even in high-dimensional or multimodal problems.

## Installation

You can install the development version of `Rtwalk` from GitHub using the `devtools` package:

```r
# If you don't have the devtools package, install it first
# install.packages("devtools")

devtools::install_github("rodrigosqrt3/Rtwalk", build_vignettes = TRUE)
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

burnin <- nrow(result$all_samples) * 0.2
samples <- result$all_samples[-(1:burnin), ]

par(mfrow = c(1, 2))
hist(samples, breaks = 50, freq = FALSE, 
     main = "Posterior Distribution", xlab = "Parameter Value")
lines(density(samples), col = "blue", lwd = 2)

plot(samples, type = 'l', col = "grey30", 
     main = "Trace Plot", xlab = "Iteration")
```

## Citation

This package is an implementation of the algorithm described in the following paper. If you use `Rtwalk` in your research, please cite the original work:

> Christen, J. A., & Fox, C. (2010). A general purpose sampling algorithm for continuous distributions (the t-walk). *Bayesian Analysis*, 5(2), 263-282. [doi:10.1214/10-BA603](https://doi.org/10.1214/10-BA603)

## License


This package is licensed under the GPL-3. See the `LICENSE` file for more details.
