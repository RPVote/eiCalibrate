#' Leave-One-Out Cross-Validation for Ecological Inference Methods
#'
#' Evaluates the out-of-sample predictive accuracy of ecological inference
#' methods by iteratively holding out each precinct, estimating
#' group-conditional candidate support from the remaining precincts, and
#' predicting the held-out precinct's aggregate vote share. The method whose
#' estimated support rates best reproduce observed precinct-level outcomes is,
#' in a predictive sense, best calibrated to the data.
#'
#' The prediction for held-out precinct \eqn{p} is:
#' \deqn{\hat{V}_p^{\text{Dem}} = \sum_{r=1}^{R} \pi_{pr} \cdot
#'   \hat{\theta}_{r,\text{Dem}}^{(-p)}}
#' where \eqn{\pi_{pr}} is the share of group \eqn{r} in precinct \eqn{p}'s
#' voting-age population, and \eqn{\hat{\theta}_{r,\text{Dem}}^{(-p)}} is
#' the estimated support for the target candidate among group \eqn{r}, fitted
#' on all precincts except \eqn{p}.
#'
#' @param data Precinct-level data frame. Must contain race count columns
#'   (`n_<race>`), race proportion columns (`pct_<race>`), candidate count
#'   columns, candidate proportion columns (`pct_<cand>`), and a totals column.
#' @param methods Character vector of methods to compare. Options:
#'   `"goodman"` (Goodman's ecological regression),
#'   `"rxc_default"` (RxC with default priors \eqn{\lambda_1 = 4, \lambda_2 = 2}),
#'   `"rxc_calibrated"` (RxC with user-specified priors),
#'   `"iter"` (iterative 2x2 EI via eiCompare).
#'   Default: `c("goodman", "rxc_default")`.
#' @param race_cols Character vector of race group base names (e.g.,
#'   `c("white", "black")`). The function expects columns named `n_<race>`
#'   and `pct_<race>` in the data.
#' @param cand_cols Character vector of candidate column names for vote counts
#'   (e.g., `c("cand_A", "cand_B")`). The function expects corresponding
#'   proportion columns named `pct_<cand>` (with `cand_` prefix replaced by
#'   `pct_`).
#' @param totals_col Name of the total-voters column. Default `"total"`.
#' @param dem_col_idx Integer. Which candidate column (1-indexed) to use as
#'   the prediction target. Default `1` (first candidate). The LOO-CV predicts
#'   this candidate's aggregate vote share in each held-out precinct.
#' @param calibrated_lambda1 Numeric. \eqn{\lambda_1} for `"rxc_calibrated"`.
#'   Required if `"rxc_calibrated"` is in `methods`.
#' @param calibrated_lambda2 Numeric. \eqn{\lambda_2} for `"rxc_calibrated"`.
#'   Defaults to `calibrated_lambda1` if not supplied.
#' @param sample Integer. MCMC sample size for RxC fits. Default `25000`.
#'   Shorter chains are appropriate for CV since only point estimates are needed.
#' @param burnin Integer. MCMC burn-in for RxC fits. Default `5000`.
#' @param thin Integer. MCMC thinning interval for RxC fits. Default `10`.
#' @param verbose Logical. Print progress messages. Default `TRUE`.
#' @param ... Additional arguments passed to [fit_ei()].
#'
#' @return A list with elements:
#' \describe{
#'   \item{precinct_errors}{Data frame with columns: `precinct` (row index),
#'     `method`, `predicted` (predicted vote share), `actual` (observed vote
#'     share), `error` (predicted - actual).}
#'   \item{summary}{Data frame with columns: `method`, `rmse`, `mae`,
#'     `n_success` (folds that produced estimates), `n_failed` (folds where
#'     fitting failed).}
#'   \item{settings}{List recording the arguments used.}
#' }
#'
#' @examples
#' # --- Fast example with Goodman ER only (no MCMC) ---
#' dat <- simulate_election(n_precincts = 30, seed = 42)
#' cv <- loo_cv(dat, methods = "goodman")
#' cv$summary
#' head(cv$precinct_errors)
#'
#' \dontrun{
#' # --- Compare Goodman ER vs RxC Default ---
#' cv2 <- loo_cv(dat, methods = c("goodman", "rxc_default"),
#'               sample = 10000, burnin = 2000, thin = 5)
#' cv2$summary
#'
#' # --- Include calibrated RxC ---
#' cv3 <- loo_cv(dat,
#'               methods = c("goodman", "rxc_default", "rxc_calibrated"),
#'               calibrated_lambda1 = 0.5, calibrated_lambda2 = 0.25,
#'               sample = 10000, burnin = 2000, thin = 5)
#' cv3$summary
#'
#' # --- Visualize results ---
#' plot_loocv(cv3, type = "summary")
#' plot_loocv(cv3, type = "precinct")
#' }
#'
#' @export
loo_cv <- function(data,
                   methods = c("goodman", "rxc_default"),
                   race_cols = c("white", "black"),
                   cand_cols = c("cand_A", "cand_B"),
                   totals_col = "total",
                   dem_col_idx = 1L,
                   calibrated_lambda1 = NULL,
                   calibrated_lambda2 = NULL,
                   sample = 25000,
                   burnin = 5000,
                   thin = 10,
                   verbose = TRUE,
                   ...) {

  # --- Input validation ---
  valid_methods <- c("goodman", "rxc_default", "rxc_calibrated", "iter")
  methods <- match.arg(methods, valid_methods, several.ok = TRUE)

  if ("rxc_calibrated" %in% methods && is.null(calibrated_lambda1)) {
    stop("calibrated_lambda1 must be provided when methods includes 'rxc_calibrated'.")
  }
  if (is.null(calibrated_lambda2)) {
    calibrated_lambda2 <- calibrated_lambda1
  }

  dem_col_idx <- as.integer(dem_col_idx)
  if (dem_col_idx < 1L || dem_col_idx > length(cand_cols)) {
    stop("dem_col_idx must be between 1 and ", length(cand_cols), ".")
  }

  # --- Derive column names ---
  race_pct_cols <- paste0("pct_", race_cols)
  cand_pct_cols <- sub("^cand_", "pct_", cand_cols)
  target_pct_col <- cand_pct_cols[dem_col_idx]

  # Validate columns exist
  required_cols <- c(race_pct_cols, cand_cols, totals_col)
  missing <- setdiff(required_cols, names(data))
  if (length(missing) > 0) {
    stop("Missing required columns in data: ", paste(missing, collapse = ", "),
         "\nExpected race proportion columns: ", paste(race_pct_cols, collapse = ", "),
         "\nExpected candidate columns: ", paste(cand_cols, collapse = ", "),
         "\nExpected totals column: ", totals_col)
  }

  n <- nrow(data)
  if (verbose) message(sprintf("LOO-CV: %d precincts, %d methods", n, length(methods)))

  # --- Main LOO loop ---
  results <- vector("list", n)

  for (i in seq_len(n)) {
    if (verbose) message(sprintf("  Fold %d/%d ...", i, n))

    train <- data[-i, , drop = FALSE]
    test  <- data[i, , drop = FALSE]
    actual <- test[[target_pct_col]]

    # Race shares for prediction weighting
    race_shares <- vapply(race_pct_cols, function(col) test[[col]],
                          numeric(1), USE.NAMES = FALSE)

    fold_rows <- vector("list", length(methods))
    names(fold_rows) <- methods

    for (meth in methods) {
      pred <- NA_real_

      if (meth == "goodman") {
        fit <- tryCatch(
          fit_ei(train, method = "goodman", race_cols = race_cols,
                 cand_cols = cand_cols, totals_col = totals_col),
          error = function(e) NULL
        )
        if (!is.null(fit)) {
          pred <- sum(race_shares * fit$estimates[, dem_col_idx])
        }
      }

      if (meth == "rxc_default") {
        fit <- tryCatch(
          fit_ei(train, method = "rxc", race_cols = race_cols,
                 cand_cols = cand_cols, totals_col = totals_col,
                 lambda1 = 4, lambda2 = 2,
                 sample = sample, burnin = burnin, thin = thin, ...),
          error = function(e) NULL
        )
        if (!is.null(fit)) {
          pred <- sum(race_shares * fit$estimates$point[, dem_col_idx])
        }
      }

      if (meth == "rxc_calibrated") {
        fit <- tryCatch(
          fit_ei(train, method = "rxc", race_cols = race_cols,
                 cand_cols = cand_cols, totals_col = totals_col,
                 lambda1 = calibrated_lambda1, lambda2 = calibrated_lambda2,
                 sample = sample, burnin = burnin, thin = thin, ...),
          error = function(e) NULL
        )
        if (!is.null(fit)) {
          pred <- sum(race_shares * fit$estimates$point[, dem_col_idx])
        }
      }

      if (meth == "iter") {
        fit <- tryCatch(
          fit_ei(train, method = "iter", race_cols = race_cols,
                 cand_cols = cand_cols, totals_col = totals_col),
          error = function(e) NULL
        )
        if (!is.null(fit)) {
          pred <- sum(race_shares * fit$estimates[, dem_col_idx])
        }
      }

      fold_rows[[meth]] <- data.frame(
        precinct  = i,
        method    = meth,
        predicted = pred,
        actual    = actual,
        error     = pred - actual,
        stringsAsFactors = FALSE
      )
    }

    results[[i]] <- do.call(rbind, fold_rows)
  }

  # --- Aggregate ---
  precinct_errors <- do.call(rbind, results)
  rownames(precinct_errors) <- NULL

  summary_df <- do.call(rbind, lapply(methods, function(m) {
    sub <- precinct_errors[precinct_errors$method == m, ]
    valid <- !is.na(sub$error)
    data.frame(
      method    = m,
      rmse      = if (any(valid)) sqrt(mean(sub$error[valid]^2)) else NA_real_,
      mae       = if (any(valid)) mean(abs(sub$error[valid])) else NA_real_,
      n_success = sum(valid),
      n_failed  = sum(!valid),
      stringsAsFactors = FALSE
    )
  }))
  rownames(summary_df) <- NULL

  if (verbose) {
    message("\nLOO-CV Summary:")
    for (i in seq_len(nrow(summary_df))) {
      message(sprintf("  %-20s  RMSE = %.4f  MAE = %.4f  (%d/%d succeeded)",
                      summary_df$method[i], summary_df$rmse[i],
                      summary_df$mae[i], summary_df$n_success[i],
                      summary_df$n_success[i] + summary_df$n_failed[i]))
    }
  }

  list(
    precinct_errors = precinct_errors,
    summary = summary_df,
    settings = list(
      methods = methods,
      race_cols = race_cols,
      cand_cols = cand_cols,
      totals_col = totals_col,
      dem_col_idx = dem_col_idx,
      n_precincts = n,
      sample = sample,
      burnin = burnin,
      thin = thin
    )
  )
}
