# Tests here are split into two groups:
#  - Pure-R helper functions (.fit_exgaussian, .pexgaussian, .split_trials,
#    .p_inhibit_r) and bundled-file lookups, which need no Stan backend and
#    always run.
#  - End-to-end ssrt_stan() fits, which are skipped unless cmdstanr or rstan
#    is installed (these are heavy Suggests, not required for the rest of the
#    package).

# --------------------------------------------------------- .fit_exgaussian

test_that(".fit_exgaussian recovers known ex-Gaussian parameters", {
  set.seed(123)
  true_mu <- 400; true_sigma <- 50; true_tau <- 100
  x <- SSRTcalc:::.rexgaussian(5000, true_mu, true_sigma, true_tau)

  fit <- SSRTcalc:::.fit_exgaussian(x)

  expect_named(fit, c("mu", "sigma", "tau"))
  expect_true(fit$mu    > 0)
  expect_true(fit$sigma > 0)
  expect_true(fit$tau   > 0)
  # method-of-moments should be within ~15% for n = 5000
  expect_equal(fit$mu,    true_mu,    tolerance = 0.15)
  expect_equal(fit$sigma, true_sigma, tolerance = 0.30)
  expect_equal(fit$tau,   true_tau,   tolerance = 0.20)
})

test_that(".fit_exgaussian errors with too few observations", {
  expect_error(SSRTcalc:::.fit_exgaussian(c(100, 200)), ">= 5")
})

# ------------------------------------------------------------ .pexgaussian

test_that(".pexgaussian is a valid, monotone CDF", {
  mu <- 400; sigma <- 50; tau <- 100
  x  <- seq(0, 1500, by = 10)
  p  <- SSRTcalc:::.pexgaussian(x, mu, sigma, tau)

  expect_true(all(p >= -1e-8 & p <= 1 + 1e-8))
  expect_true(all(diff(p) >= -1e-8))            # monotone non-decreasing
  expect_lt(SSRTcalc:::.pexgaussian(-500, mu, sigma, tau), 0.01)
  expect_gt(SSRTcalc:::.pexgaussian(2000, mu, sigma, tau), 0.99)
})

# -------------------------------------------------------------- .split_trials

test_that(".split_trials partitions go / signal-respond / inhibited correctly", {
  d <- adaptive[adaptive$SubjID == 1, ]
  s <- SSRTcalc:::.split_trials(d, "vol", "RT_exp", "soa")

  expect_named(s, c("go", "sr", "inh"))
  expect_equal(nrow(s$go) + nrow(s$sr) + nrow(s$inh),
               sum(d$vol == 0) + sum(d$vol == 1))

  # All go-trial RTs are valid
  expect_true(all(!is.na(s$go$RT_exp) & s$go$RT_exp > 0))
  # Signal-respond trials have an observed RT; inhibited trials do not
  expect_true(all(!is.na(s$sr$RT_exp)))
  expect_true(all(is.na(s$inh$RT_exp)))
})

# -------------------------------------------------------------- .p_inhibit_r

test_that(".p_inhibit_r returns a probability that decreases with SSD", {
  mu_go <- 500; sigma_go <- 80; tau_go <- 100
  mu_sp <- 200; sigma_sp <- 30; tau_sp <- 50

  p_short <- SSRTcalc:::.p_inhibit_r(  50, mu_go, sigma_go, tau_go,
                                            mu_sp, sigma_sp, tau_sp)
  p_long  <- SSRTcalc:::.p_inhibit_r( 600, mu_go, sigma_go, tau_go,
                                            mu_sp, sigma_sp, tau_sp)

  expect_true(p_short >= 0 && p_short <= 1)
  expect_true(p_long  >= 0 && p_long  <= 1)
  # Longer SSD -> the go process has more time to win -> lower P(inhibit)
  expect_true(p_long < p_short)
})

# ------------------------------------------------------------ bundled .stan

test_that("bundled Stan model files are present and non-empty", {
  for (f in c("ssrt_single.stan", "ssrt_hier.stan")) {
    path <- SSRTcalc:::.stan_model_path(f)
    expect_true(file.exists(path))
    expect_true(file.size(path) > 0)
    expect_true(endsWith(path, f))
  }
})

test_that(".stan_model_path errors informatively for an unknown file", {
  expect_error(SSRTcalc:::.stan_model_path("does_not_exist.stan"), "not found")
})

# ------------------------------------------------------------ backend detect

test_that(".detect_stan_backend errors with install instructions if neither backend is present", {
  has_cmdstanr <- requireNamespace("cmdstanr", quietly = TRUE)
  has_rstan    <- requireNamespace("rstan",    quietly = TRUE)

  if (!has_cmdstanr && !has_rstan) {
    expect_error(SSRTcalc:::.detect_stan_backend(), "No Stan backend found")
  } else {
    expect_type(SSRTcalc:::.detect_stan_backend(), "character")
  }
})

# -------------------------------------------------------- end-to-end fitting

test_that("ssrt_stan() fits a single-subject model end to end (slow, needs a backend)", {
  skip_if_not_installed("cmdstanr")
  has_cmdstan <- requireNamespace("cmdstanr", quietly = TRUE) &&
    !is.null(tryCatch(cmdstanr::cmdstan_path(), error = function(e) NULL))
  has_rstan <- requireNamespace("rstan", quietly = TRUE)
  skip_if_not(has_cmdstan || has_rstan, "no usable Stan backend installed")

  d <- adaptive[adaptive$SubjID == 1, ]

  fit <- ssrt_stan(d, chains = 1, iter = 400, warmup = 200,
                   cores = 1, n_quad = 40, seed = 1)

  expect_s3_class(fit, "ssrt_stan")
  expect_equal(fit$model_type, "single")
  expect_true("mean_ssrt" %in% fit$summary$variable)

  ssrt_row <- fit$summary[fit$summary$variable == "mean_ssrt", ]
  expect_true(ssrt_row$mean > 0)
})
