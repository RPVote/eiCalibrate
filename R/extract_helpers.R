#' Extract Row-Conditional Estimates from an eiPack RxC Fit
#'
#' Aggregates posterior draws of cell counts across precincts and normalizes by
#' row totals to produce group-conditional candidate support rates with
#' credible intervals.
#'
#' The exact MCMC column naming in `ei.MD.bayes` output can vary across
#' versions of eiPack. If extraction fails, inspect
#' `colnames(as.matrix(fit$draws$Cell.counts))` and report the format so the
#' regex pattern can be adjusted.
#'
#' @param fit An object returned by [eiPack::ei.MD.bayes()] with
#'   `ret.mcmc = TRUE`.
#' @param row_names Character vector of row (demographic group) names matching
#'   those used in the formula's RHS `cbind()`.
#' @param col_names Character vector of column (candidate) names matching the
#'   formula's LHS `cbind()`.
#'
#' @return A list with `point` (matrix of posterior means), `lower` and `upper`
#'   (2.5\% and 97.5\% credible bounds), and `draws` (a 3D array of normalized
#'   posterior draws: `[draw, row, col]`).
#'
#' @examples
#' \dontrun{
#' # Typically called internally by fit_ei(), but can be used directly:
#' dat <- simulate_election(n_precincts = 40, seed = 1)
#'
#' # Fit the RxC model directly via eiPack
#' fml <- cbind(cand_A, cand_B) ~ cbind(n_white, n_black)
#' raw_fit <- eiPack::ei.MD.bayes(fml,
#'   data = dat, total = "total",
#'   lambda1 = 4, lambda2 = 2,
#'   sample = 20000, burnin = 5000,
#'   thin = 10, ret.mcmc = TRUE
#' )
#'
#' # Extract and normalize estimates
#' est <- extract_rxc_estimates(raw_fit,
#'   row_names = c("n_white", "n_black"),
#'   col_names = c("cand_A", "cand_B")
#' )
#' est$point # posterior means (rows = groups, cols = candidates)
#' est$lower # 2.5% credible bounds
#' est$upper # 97.5% credible bounds
#' dim(est$draws) # [n_draws, 2, 2]
#'
#' # Compare point estimates to truth
#' attr(dat, "true_support")
#' }
#'
#' @export
extract_rxc_estimates <- function(fit, row_names, col_names) {
  if (is.null(fit$draws) || is.null(fit$draws$Cell.counts)) {
    stop(
      "Fit object does not contain $draws$Cell.counts. ",
      "Refit with ret.mcmc = TRUE."
    )
  }

  draws <- as.matrix(fit$draws$Cell.counts)
  cn <- colnames(draws)

  n_draws <- nrow(draws)
  n_rows <- length(row_names)
  n_cols <- length(col_names)

  agg <- array(NA_real_,
    dim = c(n_draws, n_rows, n_cols),
    dimnames = list(NULL, row_names, col_names)
  )

  for (r in seq_along(row_names)) {
    for (c in seq_along(col_names)) {
      # eiPack typically names columns "ccount.<row>.<col>.<precinct>"
      # We use a flexible regex to catch variants.
      pattern <- paste0(
        "(ccount|cell)\\.?",
        row_names[r], "\\.?", col_names[c], "\\.?"
      )
      matched <- grep(pattern, cn, value = TRUE)
      if (length(matched) == 0) {
        # Fallback: try just the row and col tokens in any order
        pattern2 <- paste0("(?=.*", row_names[r], ")(?=.*", col_names[c], ")")
        matched <- grep(pattern2, cn, value = TRUE, perl = TRUE)
      }
      if (length(matched) == 0) {
        warning(
          "No MCMC columns matched for row '", row_names[r],
          "' x col '", col_names[c],
          "'. Inspect colnames(as.matrix(fit$draws$Cell.counts))."
        )
        next
      }
      agg[, r, c] <- rowSums(draws[, matched, drop = FALSE])
    }
  }

  # Normalize each row to sum to 1 per draw
  for (d in seq_len(n_draws)) {
    for (r in seq_len(n_rows)) {
      tot <- sum(agg[d, r, ])
      if (!is.na(tot) && tot > 0) agg[d, r, ] <- agg[d, r, ] / tot
    }
  }

  point <- apply(agg, c(2, 3), mean, na.rm = TRUE)
  lower <- apply(agg, c(2, 3), stats::quantile, 0.025, na.rm = TRUE)
  upper <- apply(agg, c(2, 3), stats::quantile, 0.975, na.rm = TRUE)

  list(point = point, lower = lower, upper = upper, draws = agg)
}


#' Coerce an `ei_iter` Output to a Group-by-Candidate Matrix
#'
#' `eiCompare::ei_iter()` returns estimates in a long-format data frame whose
#' exact column names vary by version. This helper extracts point estimates
#' and reshapes them into a matrix matching the format used by [fit_ei()].
#'
#' @param iter_fit The object returned by [eiCompare::ei_iter()].
#' @param row_names Character vector of demographic group names (without the
#'   `pct_` prefix).
#' @param col_names Character vector of candidate names (e.g., `cand_A`).
#'
#' @return A numeric matrix with rows = groups, columns = candidates.
#'
#' @examples
#' \dontrun{
#' # Typically called internally by fit_ei(), but can be used directly:
#' dat <- simulate_election(n_precincts = 40, seed = 1)
#'
#' # Fit iterative 2x2 EI via eiCompare
#' iter_fit <- eiCompare::ei_iter(
#'   data       = dat,
#'   cand_cols  = c("pct_A", "pct_B"),
#'   race_cols  = c("pct_white", "pct_black"),
#'   totals_col = "total"
#' )
#'
#' # Reshape to a standard matrix
#' mat <- iter_to_matrix(iter_fit,
#'   row_names = c("white", "black"),
#'   col_names = c("cand_A", "cand_B")
#' )
#' mat
#'
#' # Compare to truth
#' attr(dat, "true_support")
#' }
#'
#' @export
iter_to_matrix <- function(iter_fit, row_names, col_names) {
  # Try common locations for the estimates table
  est_df <- NULL
  for (slot in c("estimates", "results", "summary", "race_cand_table")) {
    if (!is.null(iter_fit[[slot]])) {
      est_df <- iter_fit[[slot]]
      break
    }
  }
  if (is.null(est_df) && is.data.frame(iter_fit)) {
    est_df <- iter_fit
  }
  if (is.null(est_df)) {
    stop(
      "Could not locate estimates table in ei_iter output. ",
      "Inspect names(iter_fit)."
    )
  }

  mat <- matrix(NA_real_,
    nrow = length(row_names), ncol = length(col_names),
    dimnames = list(row_names, col_names)
  )

  # Find candidate columns and race rows; ei_iter typically uses race in row
  # names and candidate as columns, OR a long format with both.
  cand_pct_cols <- sub("^cand_", "pct_", col_names)
  race_pct_cols <- paste0("pct_", row_names)

  # Case 1: wide format with race rows and candidate columns
  if (all(cand_pct_cols %in% colnames(est_df)) ||
    all(col_names %in% colnames(est_df))) {
    target_cols <- if (all(cand_pct_cols %in% colnames(est_df))) {
      cand_pct_cols
    } else {
      col_names
    }
    # Find a column identifying race
    race_col <- intersect(
      c("race", "group", "Race", "Group"),
      colnames(est_df)
    )[1]
    if (!is.na(race_col)) {
      for (r in seq_along(row_names)) {
        idx <- which(est_df[[race_col]] %in%
          c(row_names[r], race_pct_cols[r]))
        if (length(idx) >= 1) {
          for (c in seq_along(col_names)) {
            mat[r, c] <- as.numeric(est_df[idx[1], target_cols[c]])
          }
        }
      }
      return(mat)
    }
  }

  # Case 2: long format with separate cand and race columns
  est_col <- intersect(c(
    "mean", "Mean", "estimate", "Estimate",
    "point", "value"
  ), colnames(est_df))[1]
  cand_col <- intersect(
    c("cand", "Cand", "candidate", "Candidate"),
    colnames(est_df)
  )[1]
  race_col2 <- intersect(
    c("race", "Race", "group", "Group"),
    colnames(est_df)
  )[1]

  if (!is.na(est_col) && !is.na(cand_col) && !is.na(race_col2)) {
    # Both cand and race columns exist alongside a mean/estimate column
    for (r in seq_along(row_names)) {
      for (cc in seq_along(col_names)) {
        # Match race: try pct_<name> and bare <name>
        race_match <- est_df[[race_col2]] %in%
          c(
            race_pct_cols[r], row_names[r],
            paste0("n_", row_names[r])
          )
        # Match candidate: try pct_<name>, cand_<name>, and bare <name>
        cand_match <- est_df[[cand_col]] %in%
          c(
            cand_pct_cols[cc], col_names[cc],
            sub("^cand_", "", col_names[cc])
          )
        idx <- which(race_match & cand_match)
        if (length(idx) >= 1) {
          mat[r, cc] <- as.numeric(est_df[idx[1], est_col])
        }
      }
    }
    if (!all(is.na(mat))) {
      return(mat)
    }
  }

  # Case 3: long format with combined label column (fallback)
  if (!is.na(est_col)) {
    label_col <- colnames(est_df)[1]
    for (r in seq_along(row_names)) {
      for (cc in seq_along(col_names)) {
        rx <- paste0(
          "(?i)", race_pct_cols[r], ".*", cand_pct_cols[cc],
          "|", row_names[r], ".*", col_names[cc]
        )
        idx <- grep(rx, est_df[[label_col]], perl = TRUE)
        if (length(idx) >= 1) {
          mat[r, cc] <- as.numeric(est_df[idx[1], est_col])
        }
      }
    }
  }

  if (all(is.na(mat))) {
    warning(
      "Could not parse ei_iter estimates table. ",
      "Returning empty matrix. Inspect the fit object directly."
    )
  }

  mat
}
