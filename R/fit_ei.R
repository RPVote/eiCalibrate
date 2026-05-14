#' Fit an Ecological Inference Model
#'
#' Unified wrapper around three ecological inference methods: RxC Bayesian
#' (via [eiPack::ei.MD.bayes()]), iterative 2x2 (via [eiCompare::ei_iter()]),
#' and Goodman's ecological regression (a simple linear regression).
#'
#' Defaults assume the column naming conventions produced by
#' [simulate_election()]: race columns `n_white`, `n_black`, `pct_white`,
#' `pct_black`, candidate columns `cand_A`, `cand_B`, `pct_A`, `pct_B`, and
#' `total`. Override via `race_cols`, `cand_cols`, `totals_col` for real data.
#'
#' @param data A precinct-level data frame.
#' @param method One of `"rxc"`, `"iter"`, `"goodman"`.
#' @param race_cols Character vector of column names for race **counts** (used
#'   by `"rxc"`) and proportions (used by `"iter"`). Provide both via
#'   `race_count_cols` and `race_pct_cols` if your data uses different names.
#'   By default, `race_cols = c("white", "black")` and the function looks for
#'   `n_white`, `n_black`, `pct_white`, `pct_black`.
#' @param cand_cols Character vector for candidate **counts** (`cand_A`,
#'   `cand_B`) and proportions (`pct_A`, `pct_B`).
#' @param totals_col Name of the total-voters column.
#' @param lambda1,lambda2 Dirichlet concentration hyperparameters for
#'   `ei.MD.bayes`. Lower values = less shrinkage. Defaults `lambda1 = 4`,
#'   `lambda2 = 2` match eiPack defaults.
#' @param sample,burnin,thin MCMC settings for RxC. Defaults are modest; raise
#'   for production runs.
#' @param verbose Logical; print progress.
#' @param ... Additional arguments passed to the underlying method.
#'
#' @return A list with elements:
#' \describe{
#'   \item{method}{The method used.}
#'   \item{estimates}{For `"rxc"`: a list with `point`, `lower`, `upper`,
#'     `draws`. For `"iter"` and `"goodman"`: a matrix of point estimates
#'     (rows = groups, cols = candidates).}
#'   \item{fit}{The raw fit object.}
#' }
#'
#' @examples
#' # --- Goodman's ecological regression (fast, no MCMC) ---
#' dat <- simulate_election(n_precincts = 40, seed = 1)
#' fit_g <- fit_ei(dat, method = "goodman")
#' fit_g$estimates
#' # Compare to embedded truth
#' attr(dat, "true_support")
#'
#' \dontrun{
#' # --- RxC Bayesian with default prior (requires eiPack) ---
#' fit_rxc <- fit_ei(dat, method = "rxc",
#'                   sample = 20000, burnin = 5000, thin = 10)
#' fit_rxc$estimates$point   # posterior means
#' fit_rxc$estimates$lower   # 2.5% credible bound
#' fit_rxc$estimates$upper   # 97.5% credible bound
#'
#' # --- RxC with calibrated (weaker) prior ---
#' fit_cal <- fit_ei(dat, method = "rxc",
#'                   lambda1 = 0.5, lambda2 = 0.5,
#'                   sample = 20000, burnin = 5000, thin = 10)
#' fit_cal$estimates$point
#'
#' # --- Iterative 2x2 EI (requires eiCompare) ---
#' fit_iter <- fit_ei(dat, method = "iter")
#' fit_iter$estimates
#'
#' # --- Score all methods against truth ---
#' truth <- attr(dat, "true_support")
#' score_estimates(fit_g$estimates, truth)
#' score_estimates(fit_rxc$estimates$point, truth)
#' score_estimates(fit_cal$estimates$point, truth)
#' }
#'
#' @export
fit_ei <- function(data,
                   method = c("rxc", "iter", "goodman"),
                   race_cols = c("white", "black"),
                   cand_cols = c("cand_A", "cand_B"),
                   totals_col = "total",
                   lambda1 = 4,
                   lambda2 = 2,
                   sample = 50000,
                   burnin = 10000,
                   thin = 25,
                   verbose = FALSE,
                   ...) {

  method <- match.arg(method)

  race_count_cols <- paste0("n_", race_cols)
  race_pct_cols   <- paste0("pct_", race_cols)
  cand_pct_cols   <- sub("^cand_", "pct_", cand_cols)

  if (method == "rxc") {
    if (!requireNamespace("eiPack", quietly = TRUE)) {
      stop("Package 'eiPack' is required for method = 'rxc'.")
    }

    # Build the cbind() formula programmatically
    lhs <- paste0("cbind(", paste(cand_cols, collapse = ", "), ")")
    rhs <- paste0("cbind(", paste(race_count_cols, collapse = ", "), ")")
    fml <- stats::as.formula(paste(lhs, "~", rhs))

    fit <- eiPack::ei.MD.bayes(
      formula = fml,
      data    = data,
      total   = totals_col,
      lambda1 = lambda1,
      lambda2 = lambda2,
      sample  = sample,
      burnin  = burnin,
      thin    = thin,
      ret.mcmc = TRUE,
      ...
    )

    est <- extract_rxc_estimates(fit,
                                 row_names = race_count_cols,
                                 col_names = cand_cols)
    # Rename to bare group/candidate names for downstream convenience
    rownames(est$point) <- race_cols
    rownames(est$lower) <- race_cols
    rownames(est$upper) <- race_cols

    return(list(method = "rxc", estimates = est, fit = fit,
                settings = list(lambda1 = lambda1, lambda2 = lambda2,
                                sample = sample, burnin = burnin,
                                thin = thin)))
  }

  if (method == "iter") {
    if (!requireNamespace("eiCompare", quietly = TRUE)) {
      stop("Package 'eiCompare' is required for method = 'iter'.")
    }

    fit <- eiCompare::ei_iter(
      data        = data,
      cand_cols   = cand_pct_cols,
      race_cols   = race_pct_cols,
      totals_col  = totals_col,
      verbose     = verbose,
      ...
    )

    est <- iter_to_matrix(fit, row_names = race_cols, col_names = cand_cols)
    return(list(method = "iter", estimates = est, fit = fit))
  }

  if (method == "goodman") {
    # Goodman's ER: regress candidate proportion on race proportions w/o intercept
    est <- matrix(NA_real_,
                  nrow = length(race_cols), ncol = length(cand_cols),
                  dimnames = list(race_cols, cand_cols))

    rhs_terms <- paste(race_pct_cols, collapse = " + ")

    fits <- list()
    for (cand in cand_cols) {
      pct_col <- sub("^cand_", "pct_", cand)
      fml <- stats::as.formula(paste(pct_col, "~", rhs_terms, "- 1"))
      m <- stats::lm(fml, data = data)
      fits[[cand]] <- m
      coefs <- stats::coef(m)
      for (r in seq_along(race_cols)) {
        est[r, cand] <- coefs[race_pct_cols[r]]
      }
    }

    est <- pmin(pmax(est, 0), 1)
    return(list(method = "goodman", estimates = est, fit = fits))
  }
}
