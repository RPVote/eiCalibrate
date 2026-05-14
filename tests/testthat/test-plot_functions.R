test_that("plot_calibration returns ggplot for 1D input", {
  skip_if_not_installed("ggplot2")

  mock <- list(
    summary = data.frame(
      lambda  = c(0.1, 0.5, 1, 2, 4),
      rmse    = c(0.12, 0.08, 0.10, 0.13, 0.18),
      mae     = c(0.10, 0.06, 0.08, 0.11, 0.15),
      max_err = c(0.20, 0.14, 0.18, 0.22, 0.30)
    ),
    best_lambda = 0.5,
    best_lambda1 = 0.5,
    best_lambda2 = 0.5,
    loss = "rmse",
    mode = "1D"
  )

  p <- plot_calibration(mock)
  expect_s3_class(p, "gg")
})

test_that("plot_calibration returns ggplot for 2D input", {
  skip_if_not_installed("ggplot2")

  grid <- expand.grid(lambda1 = c(0.1, 0.5, 2),
                       lambda2 = c(0.25, 1, 4))
  grid$rmse <- runif(9, 0.05, 0.25)
  grid$mae  <- grid$rmse * 0.8
  grid$max_err <- grid$rmse * 1.5

  mock <- list(
    summary = grid,
    best_lambda1 = 0.5,
    best_lambda2 = 0.25,
    loss = "rmse",
    mode = "2D"
  )

  p <- plot_calibration(mock)
  expect_s3_class(p, "gg")
})

test_that("plot_sensitivity returns ggplot for 1D input", {
  skip_if_not_installed("ggplot2")

  arr <- array(runif(8), dim = c(2, 2, 2),
               dimnames = list(c("white", "black"),
                               c("cand_A", "cand_B"),
                               c("lambda_0.5", "lambda_2")))

  mock <- list(
    point_estimates = arr,
    lambda_grid = c(0.5, 2),
    mode = "1D"
  )

  # Test string format
  p1 <- plot_sensitivity(mock, cell = "black_cand_A")
  expect_s3_class(p1, "gg")

  # Test vector format
  p2 <- plot_sensitivity(mock, cell = c("white", "cand_B"))
  expect_s3_class(p2, "gg")
})

test_that("plot_sensitivity errors on invalid cell", {
  skip_if_not_installed("ggplot2")

  arr <- array(runif(4), dim = c(2, 2, 1),
               dimnames = list(c("white", "black"),
                               c("cand_A", "cand_B"),
                               "lambda_1"))

  mock <- list(point_estimates = arr, lambda_grid = 1, mode = "1D")
  expect_error(plot_sensitivity(mock, cell = "foo_bar"), "not found")
})

test_that("plot_comparison returns ggplot", {
  skip_if_not_installed("ggplot2")

  comp <- data.frame(
    method = c("RxC (default)", "Goodman ER", "Truth"),
    white_cand_A = c(0.65, 0.78, 0.80),
    black_cand_A = c(0.30, 0.18, 0.15),
    white_cand_B = c(0.35, 0.22, 0.20),
    black_cand_B = c(0.70, 0.82, 0.85),
    stringsAsFactors = FALSE
  )

  p <- plot_comparison(comp)
  expect_s3_class(p, "gg")
})

test_that("plot_loocv returns ggplot for summary type", {
  skip_if_not_installed("ggplot2")

  mock <- list(
    precinct_errors = data.frame(
      precinct = rep(1:5, 2),
      method = rep(c("goodman", "rxc_default"), each = 5),
      predicted = runif(10, 0.3, 0.7),
      actual = runif(10, 0.3, 0.7),
      error = rnorm(10, 0, 0.05),
      stringsAsFactors = FALSE
    ),
    summary = data.frame(
      method = c("goodman", "rxc_default"),
      rmse = c(0.05, 0.08),
      mae = c(0.04, 0.06),
      n_success = c(5, 5),
      n_failed = c(0, 0),
      stringsAsFactors = FALSE
    )
  )

  p1 <- plot_loocv(mock, type = "summary")
  expect_s3_class(p1, "gg")

  p2 <- plot_loocv(mock, type = "precinct")
  expect_s3_class(p2, "gg")
})
