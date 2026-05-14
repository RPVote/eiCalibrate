## eiCalibrate Worked Examples
## ============================================================================
##
## Two scenarios demonstrating the calibration workflow:
##
##   Part A (Sections 1-5): Standard calibration with well-identified data
##   Part B (Sections 6-9): Calibration rescuing poorly-identified data
##                          (homogeneous precincts → ~50-50 estimates)
##
## Run line-by-line, or source the whole thing. Requires eiPack, eiCompare,
## and coda installed.

library(eiCalibrate)


## ============================================================================
## PART A: STANDARD CALIBRATION WORKFLOW
## ============================================================================

## ---- 1. Calibration election with known ground truth -----------------------

calib_truth <- matrix(c(0.80, 0.20,
                        0.15, 0.85),
                      nrow = 2, byrow = TRUE,
                      dimnames = list(c("white", "black"),
                                      c("cand_A", "cand_B")))

calib_data <- simulate_election(n_precincts = 40,
                                true_support = calib_truth,
                                noise = 0.05,
                                seed = 42)

head(calib_data)

## ---- 2. Sweep priors against the calibration benchmark ---------------------

calib_res <- calibrate_rxc(
  calib_data    = calib_data,
  calib_truth   = calib_truth,
  lambda_grid   = c(0.1, 0.25, 0.5, 1, 2, 4),
  loss          = "rmse",
  sample        = 20000,
  burnin        = 5000,
  thin          = 10
)

print(calib_res$summary)
cat("Best lambda by RMSE:", calib_res$best_lambda, "\n")

## ---- 3. Target election (different election, same jurisdiction) ------------

## In real work the truth is unknown; here we set it to evaluate the method.
target_truth <- matrix(c(0.75, 0.25,
                         0.20, 0.80),
                       nrow = 2, byrow = TRUE,
                       dimnames = list(c("white", "black"),
                                       c("cand_A", "cand_B")))

target_data <- simulate_election(n_precincts = 35,
                                 true_support = target_truth,
                                 noise = 0.05,
                                 seed = 99)

## ---- 4. Apply calibrated prior to target, compare to alternatives ----------

comparison <- compare_methods(
  data              = target_data,
  calibrated_lambda = calib_res$best_lambda,
  truth             = target_truth,
  sample            = 50000,
  burnin            = 10000,
  thin              = 25
)

print(comparison)

## ---- 5. Optional: convergence diagnostics on the target fit ----------------

target_fit_calib <- fit_ei(target_data, method = "rxc",
                           lambda1 = calib_res$best_lambda,
                           lambda2 = calib_res$best_lambda,
                           sample = 50000, burnin = 10000, thin = 25)

## Posterior draws array: dim [draw, group, candidate]
str(target_fit_calib$estimates$draws)

## Credible intervals
target_fit_calib$estimates$lower
target_fit_calib$estimates$upper


## ============================================================================
## PART B: CALIBRATION RESCUING POORLY-IDENTIFIED DATA
## ============================================================================
##
## THE PROBLEM:
## When precincts are racially homogeneous (similar racial composition),
## ecological inference has little variation to distinguish group-specific
## voting patterns. With the default Dirichlet prior (lambda = 4), the
## prior dominates the posterior and shrinks estimates toward the pooled
## mean (~50-50), even when the true pattern is strongly polarized (80-20).
##
## This is a real and common problem in VRA analysis: jurisdictions with
## limited residential segregation produce data that is hard to decompose,
## and the eiPack default prior can mask genuine racially polarized voting.
##
## THE SOLUTION:
## Using a separate calibration election where ground truth is known (from
## exit polls, validated surveys, or homogeneous-precinct analysis), we
## calibrate the prior strength. A weaker prior lets the (faint) data signal
## through, recovering estimates closer to the true polarization pattern.

## ---- 6. The identification problem: homogeneous precincts ------------------

## True voting pattern: strongly polarized (80-20 / 20-80)
hard_truth <- matrix(c(0.80, 0.20,
                        0.20, 0.80),
                      nrow = 2, byrow = TRUE,
                      dimnames = list(c("white", "black"),
                                      c("cand_A", "cand_B")))

## Simulate a HARD dataset: all precincts have ~50/50 racial composition.
## composition_shape = c(40, 40) gives Beta(40,40), with SD ~ 0.055,
## so nearly all precincts are between 40% and 60% white.
## Also: fewer precincts (20) and more noise (0.12) to compound the problem.
hard_data <- simulate_election(
  n_precincts      = 20,
  true_support     = hard_truth,
  precinct_size_mean = 600,
  noise            = 0.12,
  composition_shape = c(40, 40),
  seed             = 314
)

## Verify the identification problem: all precincts have similar composition
cat("Precinct racial composition (pct_white):\n")
print(summary(hard_data$pct_white))
cat("Range:", round(range(hard_data$pct_white), 3), "\n\n")

## Because both groups are ~50% of every precinct and voting is 80-20/20-80,
## overall precinct vote shares are all near 50-50:
cat("Precinct vote shares (pct_A):\n")
print(summary(hard_data$pct_A))

## ---- 7. Default EI produces ~50-50 on the hard data ------------------------

## Fit with the eiPack default prior (lambda = 4): strong shrinkage
fit_default <- fit_ei(hard_data, method = "rxc",
                      lambda1 = 4, lambda2 = 2,
                      sample = 30000, burnin = 10000, thin = 10)

cat("\n--- DEFAULT PRIOR (lambda1=4, lambda2=2) ---\n")
cat("Estimates:\n")
print(round(fit_default$estimates$point, 3))
cat("\nTrue values:\n")
print(hard_truth)
cat("\nRMSE:", round(score_estimates(fit_default$estimates$point, hard_truth)$rmse, 4), "\n")

## The default prior pulls estimates toward ~0.50 for all cells.
## The model cannot distinguish white from black voting patterns because
## precincts lack variation in composition --- the prior dominates.

## Also check Goodman ER on this hard data (it also struggles)
fit_goodman_hard <- fit_ei(hard_data, method = "goodman")
cat("\n--- GOODMAN ER ---\n")
cat("Estimates:\n")
print(round(fit_goodman_hard$estimates, 3))
cat("RMSE:", round(score_estimates(fit_goodman_hard$estimates, hard_truth)$rmse, 4), "\n")

## ---- 8. Calibration against a benchmark election ---------------------------

## In practice, the calibration election comes from a SEPARATE election in
## the same jurisdiction where ground truth is available (e.g., from exit
## polls, validated voter files, or homogeneous-precinct analysis).
##
## For this example, we simulate a calibration election with the SAME
## difficult data conditions (homogeneous precincts) but with known truth.
## The key assumption: if the calibration and target elections share similar
## data characteristics, the optimal prior strength transfers.

calib_hard_truth <- matrix(c(0.80, 0.20,
                              0.20, 0.80),
                            nrow = 2, byrow = TRUE,
                            dimnames = list(c("white", "black"),
                                            c("cand_A", "cand_B")))

calib_hard_data <- simulate_election(
  n_precincts      = 20,
  true_support     = calib_hard_truth,
  precinct_size_mean = 600,
  noise            = 0.12,
  composition_shape = c(40, 40),
  seed             = 271
)

## Sweep lambda: test a range from very weak (0.05) to the eiPack default (4)
calib_hard_res <- calibrate_rxc(
  calib_data  = calib_hard_data,
  calib_truth = calib_hard_truth,
  lambda_grid = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4),
  loss        = "rmse",
  sample      = 30000,
  burnin      = 10000,
  thin        = 10
)

cat("\n--- CALIBRATION RESULTS (hard data) ---\n")
print(calib_hard_res$summary)
cat("\nBest lambda:", calib_hard_res$best_lambda, "\n")
cat("Best fit estimates:\n")
print(round(calib_hard_res$best_fit$estimates$point, 3))

## Expect: smaller lambdas outperform larger ones because they allow
## the weak data signal through rather than drowning it in prior shrinkage.

## ---- 9. Apply calibrated prior to the target and compare -------------------

## Re-fit the target (hard) data with the calibrated prior
fit_calibrated <- fit_ei(hard_data, method = "rxc",
                         lambda1 = calib_hard_res$best_lambda,
                         lambda2 = calib_hard_res$best_lambda,
                         sample = 30000, burnin = 10000, thin = 10)

cat("\n====================================================================\n")
cat("COMPARISON: DEFAULT vs CALIBRATED on hard (homogeneous-precinct) data\n")
cat("====================================================================\n\n")

cat("True values:\n")
print(hard_truth)

cat("\nDefault prior (lambda=4, lambda=2) estimates:\n")
print(round(fit_default$estimates$point, 3))
cat("  RMSE:", round(score_estimates(fit_default$estimates$point, hard_truth)$rmse, 4), "\n")
cat("  MAE: ", round(score_estimates(fit_default$estimates$point, hard_truth)$mae, 4), "\n")

cat("\nCalibrated prior (lambda=", calib_hard_res$best_lambda, ") estimates:\n")
print(round(fit_calibrated$estimates$point, 3))
cat("  RMSE:", round(score_estimates(fit_calibrated$estimates$point, hard_truth)$rmse, 4), "\n")
cat("  MAE: ", round(score_estimates(fit_calibrated$estimates$point, hard_truth)$mae, 4), "\n")

cat("\nGoodman ER estimates:\n")
print(round(fit_goodman_hard$estimates, 3))
cat("  RMSE:", round(score_estimates(fit_goodman_hard$estimates, hard_truth)$rmse, 4), "\n")

cat("\nCredible intervals (calibrated fit):\n")
cat("  Lower:\n")
print(round(fit_calibrated$estimates$lower, 3))
cat("  Upper:\n")
print(round(fit_calibrated$estimates$upper, 3))

## ---- 10. Sensitivity sweep to visualize the prior's influence --------------

## This shows how dramatically prior strength affects estimates when
## precincts are homogeneous --- the hallmark of an identification problem.

sweep_hard <- sensitivity_sweep(
  hard_data,
  lambda_grid = c(0.05, 0.1, 0.25, 0.5, 1, 2, 4),
  sample = 30000, burnin = 10000, thin = 10
)

cat("\n--- SENSITIVITY SWEEP: white_cand_A estimate by lambda ---\n")
cat("(True value = 0.80)\n\n")
for (i in seq_along(sweep_hard$lambda_grid)) {
  lam <- sweep_hard$lambda_grid[i]
  est <- sweep_hard$point_estimates["white", "cand_A", i]
  cat(sprintf("  lambda = %5.2f  -->  white_cand_A = %.3f\n", lam, est))
}

cat("\nNotice: with large lambda, the estimate is pulled toward ~0.50.\n")
cat("With small lambda, the estimate moves closer to the true 0.80.\n")
cat("This pattern is exactly why calibration helps.\n")

## ---- Summary ---------------------------------------------------------------
##
## Key takeaways:
##
## 1. When precincts have similar racial composition (low variation in
##    pct_white/pct_black), ecological inference faces an identification
##    problem. The data alone cannot distinguish group-specific voting
##    patterns because all precincts look the same.
##
## 2. The default eiPack prior (lambda1=4, lambda2=2) is designed for
##    well-identified data. On hard data, it dominates the posterior and
##    shrinks estimates toward the pooled mean (~50-50), masking genuine
##    polarization.
##
## 3. Calibrating against a benchmark election with known truth identifies
##    a weaker prior that reduces this shrinkage. The calibrated prior lets
##    the (faint) data signal through, producing estimates closer to truth.
##
## 4. Calibration is NOT a magic fix --- with truly uninformative data,
##    credible intervals will remain wide. But it systematically reduces
##    bias from prior over-shrinkage, which is the dominant error source
##    in low-variation settings.
##
## 5. Always report:
##    (a) the calibration election and source of ground truth
##    (b) the lambda grid and loss function used
##    (c) the sensitivity of results to prior choice
##    (d) credible intervals, not just point estimates
