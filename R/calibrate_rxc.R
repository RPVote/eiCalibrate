#' Calibrate the RxC Dirichlet Prior Against Known Benchmarks
#'
#' Fits [eiPack::ei.MD.bayes()] (via [fit_ei()]) across a grid of Dirichlet
#' concentration values on a *calibration* dataset where the true
#' group-by-candidate support is known. Returns the prior setting that
#' minimizes the chosen loss against ground truth. The calibrated values are
#' then meant to be applied (via [fit_ei()]) to a *separate target election*
#' from the same jurisdiction/context.
#'
#' By default, the function performs a **2D joint calibration** over separate
#' grids for `lambda1` (candidate support concentration) and `lambda2`
#' (row-marginal concentration). If only `lambda_grid` is supplied, both
#' parameters are tied to the same value (1D sweep) for backward
#' compatibility.
#'
#' This is empirical Bayes calibration. It is only defensible when (a) the
#' calibration ground truth is genuinely credible (exit polls, validated
#' surveys, homogeneous precinct analysis), and (b) the calibration and target
#' elections share comparable polarization dynamics. **Do not** tune the prior
#' to ground truth from the *same* election you intend to analyze - that is
#' circular.
#'
#' @param calib_data Precinct-level data frame for the calibration election.
#' @param calib_truth Numeric matrix of known group-conditional candidate
#'   support (rows sum to 1).
#' @param lambda1_grid Numeric vector of `lambda1` values to evaluate. If
#'   `NULL` (default), falls back to `lambda_grid` (1D mode).
#' @param lambda2_grid Numeric vector of `lambda2` values to evaluate. If
#'   `NULL` (default), falls back to `lambda_grid` (1D mode).
#' @param lambda_grid Numeric vector for backward-compatible 1D mode where
#'   `lambda1 = lambda2` at each grid point. Ignored if both `lambda1_grid`
#'   and `lambda2_grid` are supplied.
#' @param n_reps Integer. Number of replications per grid point (to smooth
#'   MCMC noise). Point estimates are averaged across reps before scoring.
#'   Default 1.
#' @param loss One of `"rmse"`, `"mae"`, `"max_err"`. Loss function for
#'   selecting the best prior.
#' @param ... Additional arguments passed to [fit_ei()] (e.g. `sample`,
#'   `burnin`, `thin`, `race_cols`, `cand_cols`, `totals_col`).
#'
#' @return A list with:
#' \describe{
#'   \item{summary}{Data frame of results across the grid. In 2D mode,
#'     columns include `lambda1`, `lambda2`, `rmse`, `mae`, `max_err`.}
#'   \item{best_lambda1}{The optimal `lambda1`.}
#'   \item{best_lambda2}{The optimal `lambda2`.}
#'   \item{best_lambda}{(1D mode only) The single optimal lambda. In 2D mode,
#'     this is `NULL`.}
#'   \item{best_fit}{The fit (or averaged estimates) at the optimal point.}
#'   \item{rmse_matrix}{(2D mode only) Matrix of RMSE values with `lambda1`
#'     on rows and `lambda2` on columns.}
#'   \item{loss}{The loss function used.}
#'   \item{mode}{`"1D"` or `"2D"` indicating which calibration was run.}
#' }
#'
#' @examples
#' \dontrun{
#' truth <- matrix(c(0.80, 0.20, 0.15, 0.85), nrow = 2, byrow = TRUE,
#'                 dimnames = list(c("white", "black"),
#'                                 c("cand_A", "cand_B")))
#' calib <- simulate_election(n_precincts = 40, true_support = truth, seed = 1)
#'
#' # --- 2D joint calibration (recommended) ---
#' calib_2d <- calibrate_rxc(calib, truth,
#'                           lambda1_grid = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4),
#'                           lambda2_grid = c(0.05, 0.25, 0.5, 1, 2, 4, 8),
#'                           n_reps = 3,
#'                           sample = 20000, burnin = 5000, thin = 10)
#' cat("Best (lambda1, lambda2):", calib_2d$best_lambda1, ",",
#'     calib_2d$best_lambda2, "\n")
#' print(calib_2d$rmse_matrix)
#'
#' # --- 1D calibration (backward compatible) ---
#' calib_1d <- calibrate_rxc(calib, truth,
#'                           lambda_grid = c(0.1, 0.25, 0.5, 1, 2, 4),
#'                           sample = 20000, burnin = 5000, thin = 10)
#' cat("Best lambda:", calib_1d$best_lambda, "\n")
#' }
#'
#' @export
calibrate_rxc <- function(calib_data, calib_truth,
                          lambda1_grid = NULL,
                          lambda2_grid = NULL,
                          lambda_grid = c(0.1, 0.25, 0.5, 1, 2, 4),
                          n_reps = 1,
                          loss = c("rmse", "mae", "max_err"),
                          ...) {

  loss <- match.arg(loss)
  n_reps <- max(as.integer(n_reps), 1L)

  # Determine 1D vs 2D mode
  use_2d <- !is.null(lambda1_grid) && !is.null(lambda2_grid)

  if (use_2d) {
    # --- 2D joint calibration ---
    grid <- expand.grid(lambda1 = lambda1_grid, lambda2 = lambda2_grid,
                        KEEP.OUT.ATTRS = FALSE)

    results <- vector("list", nrow(grid))
    for (g in seq_len(nrow(grid))) {
      l1 <- grid$lambda1[g]
      l2 <- grid$lambda2[g]
      message(sprintf("Calibrating lambda1 = %g, lambda2 = %g (config %d/%d)",
                      l1, l2, g, nrow(grid)))

      rep_estimates <- vector("list", n_reps)
      for (rep in seq_len(n_reps)) {
        fit <- fit_ei(calib_data, method = "rxc",
                      lambda1 = l1, lambda2 = l2, ...)
        rep_estimates[[rep]] <- fit$estimates$point
      }

      # Average across replications
      avg_est <- Reduce("+", rep_estimates) / n_reps
      score <- score_estimates(avg_est, calib_truth)
      results[[g]] <- list(lambda1 = l1, lambda2 = l2,
                           avg_estimates = avg_est, score = score)
    }

    summary_df <- data.frame(
      lambda1 = grid$lambda1,
      lambda2 = grid$lambda2,
      rmse    = vapply(results, function(x) x$score$rmse,    numeric(1)),
      mae     = vapply(results, function(x) x$score$mae,     numeric(1)),
      max_err = vapply(results, function(x) x$score$max_err, numeric(1))
    )

    best_idx <- which.min(summary_df[[loss]])
    best_l1  <- summary_df$lambda1[best_idx]
    best_l2  <- summary_df$lambda2[best_idx]

    # Build RMSE matrix (lambda1 on rows, lambda2 on columns)
    rmse_mat <- matrix(summary_df$rmse, nrow = length(lambda1_grid),
                       ncol = length(lambda2_grid), byrow = FALSE)
    # Fix: reshape properly from the expand.grid ordering
    rmse_mat <- matrix(NA_real_, nrow = length(lambda1_grid),
                       ncol = length(lambda2_grid),
                       dimnames = list(as.character(lambda1_grid),
                                       as.character(lambda2_grid)))
    for (i in seq_len(nrow(summary_df))) {
      ri <- match(summary_df$lambda1[i], lambda1_grid)
      ci <- match(summary_df$lambda2[i], lambda2_grid)
      rmse_mat[ri, ci] <- summary_df$rmse[i]
    }

    list(
      summary      = summary_df,
      best_lambda1 = best_l1,
      best_lambda2 = best_l2,
      best_lambda  = NULL,
      best_fit     = results[[best_idx]]$avg_estimates,
      rmse_matrix  = rmse_mat,
      all_results  = results,
      loss         = loss,
      mode         = "2D"
    )

  } else {
    # --- 1D calibration (backward compatible) ---
    results <- lapply(lambda_grid, function(lam) {
      message("Calibrating with lambda = ", lam)

      rep_estimates <- vector("list", n_reps)
      rep_fits <- vector("list", n_reps)
      for (rep in seq_len(n_reps)) {
        fit <- fit_ei(calib_data, method = "rxc",
                      lambda1 = lam, lambda2 = lam, ...)
        rep_estimates[[rep]] <- fit$estimates$point
        rep_fits[[rep]] <- fit
      }

      avg_est <- Reduce("+", rep_estimates) / n_reps
      score <- score_estimates(avg_est, calib_truth)
      list(lambda = lam, fit = rep_fits[[1]], avg_estimates = avg_est,
           score = score)
    })

    summary_df <- data.frame(
      lambda  = lambda_grid,
      rmse    = vapply(results, function(x) x$score$rmse,    numeric(1)),
      mae     = vapply(results, function(x) x$score$mae,     numeric(1)),
      max_err = vapply(results, function(x) x$score$max_err, numeric(1))
    )

    best_idx    <- which.min(summary_df[[loss]])
    best_lambda <- summary_df$lambda[best_idx]
    best_fit    <- results[[best_idx]]$fit

    list(
      summary      = summary_df,
      best_lambda  = best_lambda,
      best_lambda1 = best_lambda,
      best_lambda2 = best_lambda,
      best_fit     = best_fit,
      rmse_matrix  = NULL,
      all_results  = results,
      loss         = loss,
      mode         = "1D"
    )
  }
}
