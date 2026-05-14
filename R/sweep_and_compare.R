#' Sweep RxC Prior Strength Without a Ground-Truth Benchmark
#'
#' Like [calibrate_rxc()] but without scoring against truth - returns point
#' estimates across a grid of prior strengths so you can assess sensitivity.
#' Use this when you have no calibration data but want to understand how much
#' the prior is driving your results.
#'
#' By default, performs a 1D sweep where `lambda1 = lambda2`. Supply separate
#' `lambda1_grid` and `lambda2_grid` for a 2D sensitivity analysis.
#'
#' @param data Precinct-level data.
#' @param lambda1_grid Numeric vector of `lambda1` values. If `NULL`,
#'   falls back to `lambda_grid` (1D mode).
#' @param lambda2_grid Numeric vector of `lambda2` values. If `NULL`,
#'   falls back to `lambda_grid` (1D mode).
#' @param lambda_grid Numeric vector for 1D mode (lambda1 = lambda2).
#' @param ... Passed to [fit_ei()].
#'
#' @return A list with `all_results` (list of fits), `point_estimates` (an
#'   array of group-by-candidate estimates), `lambda_grid` (1D) or
#'   `lambda1_grid`/`lambda2_grid` (2D), and `mode`.
#'
#' @examples
#' \dontrun{
#' dat <- simulate_election(n_precincts = 40, seed = 1)
#'
#' # 1D sweep (lambda1 = lambda2)
#' sweep_1d <- sensitivity_sweep(dat,
#'                               lambda_grid = c(0.1, 0.5, 1, 2, 4),
#'                               sample = 20000, burnin = 5000, thin = 10)
#' sweep_1d$point_estimates
#'
#' # 2D sweep (independent lambda1 and lambda2)
#' sweep_2d <- sensitivity_sweep(dat,
#'                               lambda1_grid = c(0.1, 0.5, 1, 4),
#'                               lambda2_grid = c(0.25, 1, 2, 8),
#'                               sample = 20000, burnin = 5000, thin = 10)
#' sweep_2d$point_estimates  # 4D array: [group, candidate, lambda1, lambda2]
#' }
#'
#' @export
sensitivity_sweep <- function(data,
                              lambda1_grid = NULL,
                              lambda2_grid = NULL,
                              lambda_grid = c(0.1, 0.25, 0.5, 1, 2, 4),
                              ...) {

  use_2d <- !is.null(lambda1_grid) && !is.null(lambda2_grid)

  if (use_2d) {
    # --- 2D sensitivity sweep ---
    grid <- expand.grid(lambda1 = lambda1_grid, lambda2 = lambda2_grid,
                        KEEP.OUT.ATTRS = FALSE)
    results <- vector("list", nrow(grid))
    for (g in seq_len(nrow(grid))) {
      l1 <- grid$lambda1[g]
      l2 <- grid$lambda2[g]
      message(sprintf("Fitting lambda1 = %g, lambda2 = %g (%d/%d)",
                      l1, l2, g, nrow(grid)))
      results[[g]] <- fit_ei(data, method = "rxc",
                             lambda1 = l1, lambda2 = l2, ...)
    }

    # Build a 4D array: [group, candidate, lambda1_idx, lambda2_idx]
    pe <- results[[1]]$estimates$point
    arr <- array(NA_real_,
                 dim = c(nrow(pe), ncol(pe),
                         length(lambda1_grid), length(lambda2_grid)),
                 dimnames = list(rownames(pe), colnames(pe),
                                 paste0("l1_", lambda1_grid),
                                 paste0("l2_", lambda2_grid)))
    for (g in seq_len(nrow(grid))) {
      ri <- match(grid$lambda1[g], lambda1_grid)
      ci <- match(grid$lambda2[g], lambda2_grid)
      arr[, , ri, ci] <- results[[g]]$estimates$point
    }

    list(all_results = results,
         point_estimates = arr,
         lambda1_grid = lambda1_grid,
         lambda2_grid = lambda2_grid,
         mode = "2D")

  } else {
    # --- 1D sensitivity sweep (backward compatible) ---
    results <- lapply(lambda_grid, function(lam) {
      message("Fitting with lambda = ", lam)
      fit_ei(data, method = "rxc",
             lambda1 = lam, lambda2 = lam, ...)
    })

    pe <- results[[1]]$estimates$point
    arr <- array(NA_real_,
                 dim = c(nrow(pe), ncol(pe), length(lambda_grid)),
                 dimnames = list(rownames(pe), colnames(pe),
                                 paste0("lambda_", lambda_grid)))
    for (i in seq_along(results)) {
      arr[, , i] <- results[[i]]$estimates$point
    }

    list(all_results = results,
         point_estimates = arr,
         lambda_grid = lambda_grid,
         mode = "1D")
  }
}


#' Compare RxC, ei_iter, and Goodman Estimates Side by Side
#'
#' Runs all three methods (and optionally a calibrated RxC) on the same data,
#' returning a comparison data frame of group-conditional candidate support
#' estimates.
#'
#' @param data Precinct-level data.
#' @param calibrated_lambda1 Optional numeric. `lambda1` for the calibrated
#'   RxC fit. If supplied alongside `calibrated_lambda2`, a calibrated fit is
#'   produced.
#' @param calibrated_lambda2 Optional numeric. `lambda2` for the calibrated
#'   RxC fit. Defaults to `calibrated_lambda1` if not supplied (backward
#'   compatible behavior).
#' @param calibrated_lambda Optional numeric. Backward-compatible shorthand:
#'   sets both `lambda1` and `lambda2` to this value. Ignored if
#'   `calibrated_lambda1` is supplied.
#' @param truth Optional truth matrix. If supplied, included as a row in the
#'   output for direct comparison (useful in simulations).
#' @param ... Passed to [fit_ei()] for the RxC fits.
#'
#' @return A data frame with one row per method and columns for each
#'   group-candidate cell.
#'
#' @examples
#' \dontrun{
#' truth <- matrix(c(0.80, 0.20, 0.15, 0.85), nrow = 2, byrow = TRUE,
#'                 dimnames = list(c("white", "black"),
#'                                 c("cand_A", "cand_B")))
#' dat <- simulate_election(n_precincts = 40, true_support = truth, seed = 1)
#'
#' # Without calibration
#' comp <- compare_methods(dat, truth = truth,
#'                         sample = 20000, burnin = 5000, thin = 10)
#' print(comp)
#'
#' # With jointly calibrated lambdas
#' comp2 <- compare_methods(dat,
#'                          calibrated_lambda1 = 0.5,
#'                          calibrated_lambda2 = 0.25,
#'                          truth = truth,
#'                          sample = 20000, burnin = 5000, thin = 10)
#' print(comp2)
#'
#' # Backward compatible: single lambda
#' comp3 <- compare_methods(dat,
#'                          calibrated_lambda = 0.5,
#'                          truth = truth,
#'                          sample = 20000, burnin = 5000, thin = 10)
#' print(comp3)
#' }
#'
#' @export
compare_methods <- function(data,
                            calibrated_lambda1 = NULL,
                            calibrated_lambda2 = NULL,
                            calibrated_lambda = NULL,
                            truth = NULL,
                            ...) {

  fit_default <- fit_ei(data, method = "rxc", lambda1 = 4, lambda2 = 2, ...)
  fit_iter    <- tryCatch(fit_ei(data, method = "iter"),
                          error = function(e) {
                            warning("ei_iter failed: ", conditionMessage(e))
                            NULL
                          })
  fit_good    <- fit_ei(data, method = "goodman")

  rows <- list()
  rows[["RxC (default prior)"]] <- fit_default$estimates$point

  # Resolve calibrated lambda values
  cal_l1 <- if (!is.null(calibrated_lambda1)) calibrated_lambda1 else calibrated_lambda
  cal_l2 <- if (!is.null(calibrated_lambda2)) calibrated_lambda2 else cal_l1

  if (!is.null(cal_l1)) {
    fit_calib <- fit_ei(data, method = "rxc",
                        lambda1 = cal_l1,
                        lambda2 = cal_l2, ...)
    label <- if (cal_l1 == cal_l2) {
      paste0("RxC (calibrated, lambda=", cal_l1, ")")
    } else {
      paste0("RxC (calibrated, l1=", cal_l1, ", l2=", cal_l2, ")")
    }
    rows[[label]] <- fit_calib$estimates$point
  }

  if (!is.null(fit_iter)) {
    rows[["ei_iter"]] <- fit_iter$estimates
  }
  rows[["Goodman ER"]] <- fit_good$estimates

  if (!is.null(truth)) {
    rows[["Truth"]] <- truth
  }

  # Flatten each matrix to a labelled vector and rbind
  flatten <- function(m) {
    v <- as.vector(m)
    names(v) <- as.vector(outer(rownames(m), colnames(m),
                                paste, sep = "_"))
    v
  }

  flat <- lapply(rows, flatten)
  out <- do.call(rbind, flat)
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  out <- cbind(method = rownames(out), out)
  rownames(out) <- NULL
  out
}
