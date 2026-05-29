#' Score Ecological Inference Estimates Against Ground Truth
#'
#' Computes RMSE, MAE, and max cell-wise error between an estimated
#' group-by-candidate matrix and a known truth matrix.
#'
#' @param estimate_matrix Numeric matrix (rows = groups, cols = candidates) of
#'   point estimates. For RxC output from [fit_ei()], pass `fit$estimates$point`.
#' @param truth_matrix Numeric matrix of the same shape, with rows summing to 1.
#'
#' @return A list with `rmse`, `mae`, `max_err`, and `cell_errors` (the signed
#'   element-wise differences).
#'
#' @examples
#' # --- Basic scoring ---
#' truth <- matrix(c(0.80, 0.20, 0.15, 0.85),
#'   nrow = 2, byrow = TRUE,
#'   dimnames = list(
#'     c("white", "black"),
#'     c("cand_A", "cand_B")
#'   )
#' )
#'
#' # Simulated estimate close to truth
#' set.seed(1)
#' est_good <- truth + matrix(rnorm(4, 0, 0.03), 2, 2)
#' score_estimates(est_good, truth)
#'
#' # Simulated estimate far from truth (50-50 across the board)
#' est_bad <- matrix(0.5, 2, 2,
#'   dimnames = list(
#'     c("white", "black"),
#'     c("cand_A", "cand_B")
#'   )
#' )
#' score_estimates(est_bad, truth)
#'
#' # --- Using with fit_ei output ---
#' dat <- simulate_election(n_precincts = 40, true_support = truth, seed = 1)
#' fit_g <- fit_ei(dat, method = "goodman")
#' score_estimates(fit_g$estimates, truth)
#'
#' @export
score_estimates <- function(estimate_matrix, truth_matrix) {
  if (!all(dim(estimate_matrix) == dim(truth_matrix))) {
    stop("estimate_matrix and truth_matrix must have the same dimensions.")
  }
  err <- estimate_matrix - truth_matrix
  list(
    rmse        = sqrt(mean(err^2, na.rm = TRUE)),
    mae         = mean(abs(err), na.rm = TRUE),
    max_err     = max(abs(err), na.rm = TRUE),
    cell_errors = err
  )
}
