# SSRTcalc

Tools to estimate stop-signal reaction time (SSRT) in R, following the
consensus guidance of [Verbruggen et al. (2019, *eLife*)](https://doi.org/10.7554/eLife.46323).

**Version 2.1.0** adds three major extensions on top of the original point
estimators:

1. **Monte Carlo methods** -- bootstrap confidence intervals, parametric
   ex-Gaussian simulation, minimum-trial-count / power analysis, and
   robustness checks under violations of the horse-race assumptions.
2. **Bayesian estimation via Stan** -- single-subject and hierarchical
   ex-Gaussian horse-race models, an optional trigger-failure parameter
   ([Matzke et al., 2013](https://doi.org/10.1037/a0030543)), posterior
   inhibition functions, and posterior predictive checks.
3. **Convenience wrappers** for running a full battery of analyses
   (`run_all_mc()`) or comparing models (`ssrt_stan_compare()`) in one call.

## Installation

```r
install.packages("devtools")
devtools::install_github("agleontyev/SSRTcalc")
#Alternative
pak::pkg_install("agleontyev/SSRTcalc")
```

The Monte Carlo functions only need base R (plus `MASS`, used automatically if available). The Bayesian functions additionally need **one** of:

```r
# Recommended
install.packages("cmdstanr",
  repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()

# Alternative
install.packages("rstan")
```

## Quick start

```r
library(SSRTcalc)
data(adaptive)
d <- adaptive[adaptive$SubjID == 1, ]

# Point estimates
integration_adaptiveSSD(d)
mean_adaptiveSSD(d)
```

## Monte Carlo extensions

```r
# Bootstrap confidence interval
b <- ssrt_boot(d, n_iter = 2000)
print(b)
plot(b)

# Parametric ex-Gaussian simulation
s <- ssrt_simulate(d, n_iter = 2000)
print(s)

# How many stop trials do you need?
p <- ssrt_power(d, trial_counts = c(10, 20, 30, 50, 100), n_iter = 500)
print(p)
plot(p)

# Sensitivity to assumption violations
r <- ssrt_robustness(d, violation = "trigger_failure", n_iter = 500)
print(r)
plot(r)

# Or run all four at once
res <- run_all_mc(d, n_iter = 1000)
```

## Bayesian estimation via Stan

```r
# Single subject
fit <- ssrt_stan(d, chains = 4, iter = 2000)
print(fit)
plot(fit)

ssrt_stan_pp_check(fit)
ssrt_stan_inhibition_fn(fit)

# With a trigger-failure parameter (Matzke et al., 2013)
fit_tf <- ssrt_stan(d, trigger_failure = TRUE, adapt_delta = 0.99)

# Compare the base and trigger-failure models
cmp <- ssrt_stan_compare(d, chains = 4, iter = 2000)
cmp$comparison

# Hierarchical model across all subjects
fit_h <- ssrt_stan(adaptive, model = "hierarchical",
                    subject_col = "SubjID", chains = 4, cores = 4)
ranef(fit_h)
plot(fit_h)
```

## Data format

All functions expect one row per trial, with these default columns
(configurable via `stop_col`, `rt_col`, `acc_col`, `ssd_col`):

| Column   | Meaning                                              |
|----------|------------------------------------------------------|
| `vol`    | `0` = go trial, `1` = stop trial                      |
| `RT_exp` | Reaction time (ms); `NA` if the stop trial was inhibited |
| `correct`| Accuracy (`1` = correct / successful inhibition)      |
| `soa`    | Stop-signal delay (ms); `NA` on go trials              |

Two datasets are bundled: `adaptive` (20 subjects x 200 trials, a
staircase/adaptive-SSD design) and `fixed` (50 subjects, ~576 trials each, a fixed-SSD motion-discrimination task). Both follow the format above and are used throughout the documentation examples.

## References

Verbruggen, F., Aron, A. R., Band, G. P. H., Beste, C., Bissett, P. G.,
Brockett, A. T., ... Boehler, C. N. (2019). A consensus guide to capturing
the ability to inhibit actions and impulses: the stop-signal task. *eLife*,
8, e46323. https://doi.org/10.7554/eLife.46323

Matzke, D., Dolan, C. V., Logan, G. D., Brown, S. D., & Wagenmakers, E.-J.
(2013). Bayesian parametric estimation of stop-signal reaction time
distributions. *Journal of Experimental Psychology: General*, 142(4),
1047-1073. https://doi.org/10.1037/a0030543

## License

GPL-3
