# eiCalibrate

Empirical Bayes calibration for RxC ecological inference. Tune the Dirichlet
prior in `eiPack::ei.MD.bayes()` against known benchmarks from a calibration
election, then apply the calibrated prior to your target election.

Includes unified wrappers for:

- `eiPack` RxC Bayesian (`ei.MD.bayes`)
- `eiCompare` iterative 2x2 (`ei_iter`)
- Goodman's ecological regression

## Why

The default Dirichlet prior in `ei.MD.bayes()` (`lambda1 = 4`, `lambda2 = 2`)
can over-shrink estimates toward the pooled mean when precinct-level data is
limited or homogeneous, producing apparently non-polarized results even when
voting is in fact polarized. Weakening the prior reduces shrinkage but, on
weakly identified data, also increases variance. The principled fix is to
**calibrate** the prior strength against ground truth from a separate
benchmark election, then apply that setting to your target election.

## Install

From a local copy of the source tree:

```r
# install.packages("devtools")
devtools::install("path/to/eiCalibrate")
```

Or from a tarball / zip:

```r
install.packages("eiCalibrate_0.1.0.tar.gz", repos = NULL, type = "source")
```

Dependencies you need installed separately:

```r
install.packages(c("eiPack", "eiCompare", "coda"))
```

## Quick start

```r
library(eiCalibrate)

# 1. Calibration election with known truth (e.g. from exit polls)
calib_truth <- matrix(c(0.80, 0.20,
                        0.15, 0.85),
                      nrow = 2, byrow = TRUE,
                      dimnames = list(c("white", "black"),
                                      c("cand_A", "cand_B")))

calib_data <- simulate_election(n_precincts = 40,
                                true_support = calib_truth,
                                seed = 42)

# 2. Sweep prior strengths against the benchmark
calib_res <- calibrate_rxc(calib_data, calib_truth,
                           lambda_grid = c(0.1, 0.25, 0.5, 1, 2, 4))
calib_res$summary
calib_res$best_lambda

# 3. Apply the calibrated prior to your target election
target_fit <- fit_ei(target_data, method = "rxc",
                     lambda1 = calib_res$best_lambda,
                     lambda2 = calib_res$best_lambda)

# 4. Compare to default RxC, ei_iter, and Goodman ER as a robustness check
compare_methods(target_data,
                calibrated_lambda = calib_res$best_lambda)
```

See `inst/examples/worked_example.R` for the full pipeline.

## Methodological caveats

This package implements a tool; it does not replace methodological judgment.

- **Calibrate on a different election than the target.** Tuning the prior to
  ground truth from the same election you intend to report is circular and
  not defensible.
- **The calibration is only as good as the benchmark.** Exit polls, validated
  voter files, and homogeneous-precinct analyses are reasonable benchmarks;
  cherry-picked subsets are not.
- **Calibrated prior settings transfer best when elections are similar** in
  time, geography, office, and political context. Crossing a 2020 federal
  race with a 2024 municipal race is risky.
- **Report the calibration procedure transparently** in any expert report or
  paper. "We selected lambda = X to minimize RMSE against [benchmark]" is
  defensible; hiding the step is not.
- **Always report alongside `ei_iter` and the method of bounds.** If they
  disagree, the data don't support a strong conclusion either way.
- **If calibration RMSE is still large at the best lambda**, the underlying
  data don't identify the quantities you want and no prior tuning fixes that.

## License

MIT
