test_that("integration_adaptiveSSD returns a plausible SSRT", {
  d <- adaptive[adaptive$SubjID == 1, ]
  ssrt <- integration_adaptiveSSD(d)

  expect_type(ssrt, "double")
  expect_length(ssrt, 1)
  expect_true(is.finite(ssrt))
  # SSRTs in healthy adults are typically 150-300 ms
  expect_gt(ssrt, 50)
  expect_lt(ssrt, 500)
})

test_that("mean_adaptiveSSD returns a plausible SSRT", {
  d <- adaptive[adaptive$SubjID == 1, ]
  ssrt <- mean_adaptiveSSD(d)

  expect_type(ssrt, "double")
  expect_true(is.finite(ssrt))
  expect_gt(ssrt, 50)
  expect_lt(ssrt, 500)
})

test_that("mean method >= integration method (typical ordering)", {
  # The mean method is known to overestimate SSRT relative to the
  # integration method under typical horse-race assumptions.
  for (id in unique(adaptive$SubjID)[1:5]) {
    d <- adaptive[adaptive$SubjID == id, ]
    expect_gte(mean_adaptiveSSD(d), integration_adaptiveSSD(d) - 1e-8)
  }
})

test_that("integration_fixedSSD runs and returns a finite value", {
  d <- adaptive[adaptive$SubjID == 1, ]
  # adaptive data has multiple SSD values; fixedSSD should average across them
  ssrt <- integration_fixedSSD(d)
  expect_type(ssrt, "double")
  expect_true(is.finite(ssrt))
})

test_that("mean_fixedSSD is identical to mean_adaptiveSSD", {
  d <- adaptive[adaptive$SubjID == 1, ]
  expect_equal(mean_fixedSSD(d), mean_adaptiveSSD(d))
})

test_that("input validation catches missing columns", {
  d <- adaptive[adaptive$SubjID == 1, ]
  expect_error(integration_adaptiveSSD(d, stop_col = "nope"), "not found")
})

test_that("input validation catches non-data.frame input", {
  expect_error(integration_adaptiveSSD(matrix(1:4, 2)), "data.frame")
})

test_that("integration_adaptiveSSD errors with too few go trials", {
  d <- adaptive[adaptive$SubjID == 1, ][1:3, ]
  expect_error(integration_adaptiveSSD(d), "Fewer than 5")
})

test_that("results are consistent across all 20 subjects", {
  ssrts <- vapply(unique(adaptive$SubjID), function(id) {
    integration_adaptiveSSD(adaptive[adaptive$SubjID == id, ])
  }, numeric(1))

  expect_length(ssrts, 20)
  expect_true(all(is.finite(ssrts)))
  expect_true(all(ssrts > 0))
})
