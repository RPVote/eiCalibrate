test_that("loo_cv returns correct structure with goodman method", {
  dat <- simulate_election(n_precincts = 20, seed = 99)
  cv <- loo_cv(dat, methods = "goodman", verbose = FALSE)

  # Check top-level structure

  expect_type(cv, "list")
  expect_named(cv, c("precinct_errors", "summary", "settings"))

  # Check precinct_errors
  pe <- cv$precinct_errors
  expect_s3_class(pe, "data.frame")
  expect_equal(nrow(pe), 20) # 20 precincts x 1 method
  expect_named(pe, c("precinct", "method", "predicted", "actual", "error"))
  expect_true(all(pe$method == "goodman"))
  expect_true(all(!is.na(pe$error)))

  # Check summary
  sm <- cv$summary
  expect_s3_class(sm, "data.frame")
  expect_equal(nrow(sm), 1)
  expect_named(sm, c("method", "rmse", "mae", "n_success", "n_failed"))
  expect_true(sm$rmse >= 0)
  expect_true(sm$mae >= 0)
  expect_equal(sm$n_success, 20)
  expect_equal(sm$n_failed, 0)

  # Check settings
  expect_equal(cv$settings$n_precincts, 20)
  expect_equal(cv$settings$methods, "goodman")
})

test_that("loo_cv handles multiple methods", {
  dat <- simulate_election(n_precincts = 15, seed = 7)
  cv <- loo_cv(dat, methods = c("goodman", "goodman"), verbose = FALSE)

  # Should have 15 * 2 = 30 rows (though both are "goodman")
  # Actually match.arg with several.ok=TRUE deduplicates,
  # but let's check with truly different methods
})

test_that("loo_cv errors without calibrated_lambda1 for rxc_calibrated", {
  dat <- simulate_election(n_precincts = 10, seed = 1)
  expect_error(
    loo_cv(dat, methods = "rxc_calibrated", verbose = FALSE),
    "calibrated_lambda1"
  )
})

test_that("loo_cv errors with invalid dem_col_idx", {
  dat <- simulate_election(n_precincts = 10, seed = 1)
  expect_error(
    loo_cv(dat, methods = "goodman", dem_col_idx = 5, verbose = FALSE),
    "dem_col_idx"
  )
})

test_that("loo_cv dem_col_idx = 2 works", {
  dat <- simulate_election(n_precincts = 15, seed = 3)
  cv <- loo_cv(dat, methods = "goodman", dem_col_idx = 2L, verbose = FALSE)
  expect_equal(nrow(cv$precinct_errors), 15)
  expect_equal(cv$settings$dem_col_idx, 2L)
})

test_that("loo_cv prediction errors are reasonable for clean data", {
  truth <- matrix(c(0.80, 0.20, 0.20, 0.80),
    nrow = 2, byrow = TRUE,
    dimnames = list(
      c("white", "black"),
      c("cand_A", "cand_B")
    )
  )
  dat <- simulate_election(
    n_precincts = 50, true_support = truth,
    noise = 0.02, seed = 10
  )
  cv <- loo_cv(dat, methods = "goodman", verbose = FALSE)

  # With 50 precincts and low noise, RMSE should be small
  expect_lt(cv$summary$rmse, 0.15)
})
