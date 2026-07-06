# Small n_iter throughout keeps these tests fast while still exercising the
# full code path (resampling, simulation, plotting).

local_pdf <- function() {
  path <- tempfile(fileext = ".pdf")
  pdf(path)
  withr_like_on_exit <- function() {
    dev.off()
    unlink(path)
  }
  withr_like_on_exit
}

# ---------------------------------------------------------------- ssrt_boot

test_that("ssrt_boot returns a sensible bootstrap distribution", {
  d <- adaptive[adaptive$SubjID == 1, ]
  b <- ssrt_boot(d, n_iter = 100, seed = 1)

  expect_s3_class(b, "ssrt_boot")
  expect_equal(b$n_iter, 100)
  expect_true(b$n_valid > 0)
  expect_length(b$distribution, b$n_valid)
  expect_true(b$ci[1] <= b$mean)
  expect_true(b$ci[2] >= b$mean)
  expect_equal(b$ssrt_obs, integration_adaptiveSSD(d))
})

test_that("ssrt_boot print and plot methods run without error", {
  d <- adaptive[adaptive$SubjID == 1, ]
  b <- ssrt_boot(d, n_iter = 50, seed = 1)

  expect_output(print(b), "Bootstrap")

  finish <- local_pdf()
  expect_silent(plot(b))
  finish()
})

test_that("ssrt_boot supports an alternative ssrt_fn", {
  d <- adaptive[adaptive$SubjID == 1, ]
  b <- ssrt_boot(d, ssrt_fn = "mean_adaptiveSSD", n_iter = 50, seed = 1)

  expect_equal(b$method, "mean_adaptiveSSD")
  expect_equal(b$ssrt_obs, mean_adaptiveSSD(d))
})

# ------------------------------------------------------------ ssrt_simulate

test_that("ssrt_simulate fits an ex-Gaussian and returns a CI", {
  d <- adaptive[adaptive$SubjID == 1, ]

  expect_output(
    s <- ssrt_simulate(d, n_iter = 50, seed = 1),
    "Ex-Gaussian fit"
  )

  expect_s3_class(s, "ssrt_simulate")
  expect_true(all(c("mu", "sigma", "tau") %in% names(s$exg_params)))
  expect_true(s$exg_params$sigma > 0)
  expect_true(s$exg_params$tau   > 0)
  expect_true(s$ci[1] <= s$ci[2])
})

test_that("ssrt_simulate reports recovery_bias only when ssrt_true is given", {
  d <- adaptive[adaptive$SubjID == 1, ]

  s1 <- suppressMessages(ssrt_simulate(d, n_iter = 30, seed = 1))
  expect_null(s1$recovery_bias)

  s2 <- suppressMessages(ssrt_simulate(d, n_iter = 30, seed = 1, ssrt_true = 200))
  expect_type(s2$recovery_bias, "double")
})

test_that("ssrt_simulate print and plot methods run without error", {
  d <- adaptive[adaptive$SubjID == 1, ]
  s <- suppressMessages(ssrt_simulate(d, n_iter = 30, seed = 1))

  expect_output(print(s), "Parametric MC")

  finish <- local_pdf()
  expect_silent(plot(s))
  finish()
})

# --------------------------------------------------------------- ssrt_power

test_that("ssrt_power produces a monotone-ish precision table", {
  d <- adaptive[adaptive$SubjID == 1, ]

  p <- suppressMessages(
    ssrt_power(d, trial_counts = c(10, 50), n_iter = 20, seed = 1)
  )

  expect_s3_class(p, "ssrt_power")
  expect_equal(nrow(p$power_table), 2)
  expect_true(all(p$power_table$se > 0))
  # More trials -> generally smaller SE (allow some MC noise with n_iter=20)
  expect_true(p$power_table$se[p$power_table$n_stop_trials == 50] <
                p$power_table$se[p$power_table$n_stop_trials == 10] * 2)
})

test_that("ssrt_power print and plot methods run without error", {
  d <- adaptive[adaptive$SubjID == 1, ]
  p <- suppressMessages(
    ssrt_power(d, trial_counts = c(10, 30), n_iter = 15, seed = 1)
  )

  expect_output(print(p), "Power Analysis")

  finish <- local_pdf()
  expect_silent(plot(p))
  finish()
})

# ---------------------------------------------------------- ssrt_robustness

test_that("ssrt_robustness sweeps trigger-failure rate and returns bias estimates", {
  d <- adaptive[adaptive$SubjID == 1, ]

  r <- suppressMessages(
    ssrt_robustness(
      d,
      violation     = "trigger_failure",
      n_iter        = 20,
      trigger_range = seq(0, 0.2, by = 0.1),
      seed          = 1
    )
  )

  expect_s3_class(r, "ssrt_robustness")
  expect_true(is.null(r$correlation))
  expect_true(is.null(r$go_shift))
  expect_false(is.null(r$trigger_failure))
  expect_equal(nrow(r$trigger_failure), 3)
  expect_true(is.finite(r$baseline_ssrt))

  # Higher trigger-failure rates should bias SSRT downward (toward more
  # negative bias) on average -- check monotonic-ish trend rather than
  # exact values given MC noise at n_iter = 20.
  bias <- r$trigger_failure$bias_from_baseline
  expect_true(bias[3] <= bias[1] + 30)  # generous tolerance for noise
})

test_that("ssrt_robustness print and plot methods run without error", {
  d <- adaptive[adaptive$SubjID == 1, ]
  r <- suppressMessages(
    ssrt_robustness(d, violation = "go_shift", n_iter = 15,
                    shift_range = seq(0, 20, by = 20), seed = 1)
  )

  expect_output(print(r), "Robustness")

  finish <- local_pdf()
  expect_silent(plot(r))
  finish()
})

test_that("ssrt_robustness validates the violation argument", {
  d <- adaptive[adaptive$SubjID == 1, ]
  expect_error(
    ssrt_robustness(d, violation = "not_a_thing", n_iter = 5),
    "should be one of"
  )
})

# ------------------------------------------------------------------ run_all_mc

test_that("run_all_mc returns all four result types", {
  d <- adaptive[adaptive$SubjID == 1, ]

  res <- suppressMessages(suppressWarnings(
    run_all_mc(d, n_iter = 20, seed = 1)
  ))

  expect_named(res, c("bootstrap", "simulation", "power", "robustness"))
  expect_s3_class(res$bootstrap,  "ssrt_boot")
  expect_s3_class(res$simulation, "ssrt_simulate")
  expect_s3_class(res$power,      "ssrt_power")
  expect_s3_class(res$robustness, "ssrt_robustness")
})
