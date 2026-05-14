#' eiCalibrate: Empirical Bayes Calibration for Ecological Inference
#'
#' Tools for tuning the Dirichlet prior in RxC ecological inference against
#' known benchmarks from a calibration election, then applying the calibrated
#' prior to a target election. Includes unified wrappers for `eiPack` RxC,
#' `eiCompare` iterative 2x2, and Goodman's ER.
#'
#' @section Main functions:
#' \describe{
#'   \item{[simulate_election()]}{Generate synthetic precinct data with known
#'     truth.}
#'   \item{[fit_ei()]}{Unified wrapper around the three EI methods.}
#'   \item{[score_estimates()]}{Compute RMSE/MAE against ground truth.}
#'   \item{[calibrate_rxc()]}{Sweep RxC priors against a calibration benchmark.}
#'   \item{[sensitivity_sweep()]}{Prior sensitivity without ground truth.}
#'   \item{[compare_methods()]}{Side-by-side method comparison.}
#'   \item{[loo_cv()]}{Leave-one-out cross-validation for comparing methods.}
#'   \item{[plot_calibration()]}{Visualize calibration RMSE surface.}
#'   \item{[plot_sensitivity()]}{Visualize sensitivity sweep results.}
#'   \item{[plot_comparison()]}{Forest/dot plot of method estimates.}
#'   \item{[plot_loocv()]}{Visualize LOO-CV results.}
#' }
#'
#' @section Methodological warning:
#' Calibrating the prior is **only** defensible when the calibration ground
#' truth is genuinely credible and comes from a separate election with
#' comparable dynamics to the target. Calibrating on the same election you
#' analyze is circular. Always report the calibration procedure transparently.
#'
#' @keywords internal
"_PACKAGE"
