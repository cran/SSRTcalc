# SSRTcalc 2.1.1

## Data

* Refined the `fixed` dataset: excluded subject 2 (no recorded responses on
  any trial) and one invalid stop trial (recorded reaction time of 0 ms),
  leaving 50 subjects and 28,799 trials. Go-trial `soa` is now `NA`, matching
  the `adaptive` dataset convention. See `?fixed`.


## New features

### Monte Carlo extensions (`R/mc_extensions.R`)

* `ssrt_boot()` -- nonparametric bootstrap confidence intervals for any
  SSRT estimator (`integration_*` or `mean_*`), with `print()` and `plot()`
  methods.
* `ssrt_simulate()` -- parametric Monte Carlo via ex-Gaussian simulation of
  the horse-race model, with optional parameter-recovery checking via
  `ssrt_true`.
* `ssrt_power()` -- minimum-trial-count / power analysis: simulates datasets
  across a range of stop-trial counts and reports SE and RMSE of the SSRT
  estimate as a function of sample size.
* `ssrt_robustness()` -- sensitivity analysis for three violations of the
  independent horse-race assumptions: go/stop process correlation, trigger
  failures, and go-RT slowing on stop trials ("stop-signal cost").
* `run_all_mc()` -- convenience wrapper running all four analyses above.

### Bayesian / Stan extensions (`R/stan_extensions.R`, `inst/stan/`)

* `ssrt_stan()` -- fits the independent horse-race model with ex-Gaussian
  go- and stop-process distributions via Hamiltonian Monte Carlo.
  - `model = "single"` (default) or `model = "hierarchical"`
    (non-centred parameterisation across subjects).
  - `trigger_failure = TRUE` adds the trigger-failure parameter of
    Matzke et al. (2013).
  - Works with either `cmdstanr` (recommended) or `rstan` as backend.
* `print()`, `summary()`, `plot()`, `coef()` methods for `ssrt_stan` objects,
  including posterior densities, trace plots, pairs plots (single-subject)
  and a forest plot of per-subject SSRT (hierarchical).
* `ranef()` -- new S3 generic (SSRTcalc previously had no dependency that
  provided this), with a method for hierarchical `ssrt_stan` fits returning
  per-subject posterior summaries.
* `ssrt_stan_pp_check()` -- posterior predictive check plots for go RT and
  SSRT.
* `ssrt_stan_inhibition_fn()` -- posterior inhibition function P(inhibit |
  SSD) with a 90% credible band.
* `ssrt_stan_compare()` -- fits the base and trigger-failure models
  side by side and tabulates posterior SSRT summaries.
* `ssrt_stan_loo()` -- LOO-CV via the `loo` package (requires `log_lik` in
  generated quantities; the function explains how to add it if absent).
* Two bundled Stan programs in `inst/stan/`: `ssrt_single.stan` and
  `ssrt_hier.stan`.

### Internals

* `R/helpers.R` -- shared internal utilities: ex-Gaussian fitting
  (`.fit_exgaussian`), sampling (`.rexgaussian`), CDF (`.pexgaussian`),
  trial splitting (`.split_trials`), numerical inhibition-probability
  quadrature (`.p_inhibit_r`), and posterior-draw extraction helpers shared
  between the MC and Stan extensions.
* `R/data.R` -- documentation for the bundled `adaptive` dataset.

## Other changes

* Added a full `testthat` (edition 3) test suite covering the core
  estimators, all four Monte Carlo functions, and the pure-R Stan helpers
  (end-to-end Stan fitting tests are skipped automatically if neither
  `cmdstanr` nor `rstan` is installed).
* `Imports` now includes `graphics` and `grDevices` (used by the new
  `plot.*` methods).
* Fixed package metadata (NAMESPACE imports, ORCID).

---

# SSRTcalc 0.3.3

* Initial CRAN release: `integration_adaptiveSSD()`, `integration_fixedSSD()`,
  `mean_adaptiveSSD()`, `mean_fixedSSD()`.
