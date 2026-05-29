# --------------------------------------------------------------------------
#  Visualization functions for eiCalibrate
#  All functions require ggplot2 (in Suggests) and return ggplot objects
# --------------------------------------------------------------------------

#' Plot Calibration Results
#'
#' Visualize the RMSE (or other loss) surface from [calibrate_rxc()].
#' In 1D mode, produces a line plot of loss versus lambda. In 2D mode,
#' produces a heatmap of the loss surface across the
#' \eqn{(\lambda_1, \lambda_2)} grid with the optimal cell highlighted.
#'
#' @param calib_result Output of [calibrate_rxc()].
#' @param loss Character. Which loss to plot: `"rmse"`, `"mae"`, or
#'   `"max_err"`. Defaults to the loss used during calibration
#'   (`calib_result$loss`).
#' @param title Optional character string for the plot title.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' truth <- matrix(c(0.80, 0.20, 0.15, 0.85),
#'   nrow = 2, byrow = TRUE,
#'   dimnames = list(
#'     c("white", "black"),
#'     c("cand_A", "cand_B")
#'   )
#' )
#' dat <- simulate_election(n_precincts = 40, true_support = truth, seed = 1)
#'
#' # 1D calibration + plot
#' cal_1d <- calibrate_rxc(dat, truth,
#'   lambda_grid = c(0.1, 0.5, 1, 2, 4),
#'   sample = 10000, burnin = 2000, thin = 5
#' )
#' plot_calibration(cal_1d)
#'
#' # 2D calibration + heatmap
#' cal_2d <- calibrate_rxc(dat, truth,
#'   lambda1_grid = c(0.1, 0.5, 1, 4),
#'   lambda2_grid = c(0.25, 1, 2, 8),
#'   sample = 10000, burnin = 2000, thin = 5
#' )
#' plot_calibration(cal_2d)
#' }
#'
#' @export
plot_calibration <- function(calib_result, loss = NULL, title = NULL) {
  check_ggplot2()

  if (is.null(calib_result$mode) || is.null(calib_result$summary)) {
    stop("Input does not appear to be output from calibrate_rxc().")
  }

  loss <- if (!is.null(loss)) loss else calib_result$loss
  df <- calib_result$summary

  if (calib_result$mode == "1D") {
    # --- 1D: line plot ---
    best_lam <- calib_result$best_lambda
    best_row <- df[df$lambda == best_lam, ]

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$lambda, y = .data[[loss]])) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::geom_vline(
        xintercept = best_lam, linetype = "dashed",
        color = "red", alpha = 0.6
      ) +
      ggplot2::geom_point(data = best_row, size = 4, color = "red") +
      ggplot2::scale_x_log10() +
      ggplot2::labs(
        x = expression(lambda ~ "(log scale)"),
        y = toupper(loss),
        title = title %||% paste0("Calibration: ", toupper(loss), " vs. Lambda"),
        subtitle = sprintf("Optimal: lambda = %g", best_lam)
      ) +
      ggplot2::theme_minimal(base_size = 12)
  } else {
    # --- 2D: heatmap ---
    best_l1 <- calib_result$best_lambda1
    best_l2 <- calib_result$best_lambda2

    df$lambda1_f <- factor(df$lambda1)
    df$lambda2_f <- factor(df$lambda2)

    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = .data$lambda2_f,
      y = .data$lambda1_f,
      fill = .data[[loss]]
    )) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::scale_fill_viridis_c(
        option = "inferno", direction = -1,
        name = toupper(loss)
      )

    # Add text values if grid is manageable
    if (nrow(df) <= 64) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.3f", .data[[loss]])),
        color = "white", size = 3
      )
    }

    # Mark optimal cell
    opt_df <- data.frame(
      lambda1_f = factor(best_l1, levels = levels(df$lambda1_f)),
      lambda2_f = factor(best_l2, levels = levels(df$lambda2_f))
    )
    p <- p +
      ggplot2::geom_point(
        data = opt_df, inherit.aes = FALSE,
        ggplot2::aes(x = .data$lambda2_f, y = .data$lambda1_f),
        shape = 4, size = 6, color = "cyan", stroke = 2
      ) +
      ggplot2::labs(
        x = expression(lambda[2]),
        y = expression(lambda[1]),
        title = title %||% paste0("2D Calibration: ", toupper(loss), " Surface"),
        subtitle = sprintf("Optimal: (lambda1, lambda2) = (%g, %g)", best_l1, best_l2)
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
  }

  p
}


#' Plot Sensitivity Sweep Results
#'
#' Visualize how a specific group-conditional support estimate varies across
#' prior strength, using output from [sensitivity_sweep()].
#'
#' @param sweep_result Output of [sensitivity_sweep()].
#' @param cell Character identifying the cell to plot. Either a single string
#'   in `"group_candidate"` format (e.g., `"white_cand_A"`) or a length-2
#'   vector `c("white", "cand_A")`.
#' @param title Optional plot title.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' dat <- simulate_election(n_precincts = 40, seed = 1)
#' sweep <- sensitivity_sweep(dat,
#'   lambda_grid = c(0.1, 0.5, 1, 2, 4),
#'   sample = 10000, burnin = 2000, thin = 5
#' )
#' plot_sensitivity(sweep, cell = "black_cand_A")
#' plot_sensitivity(sweep, cell = c("white", "cand_B"))
#' }
#'
#' @export
plot_sensitivity <- function(sweep_result, cell, title = NULL) {
  check_ggplot2()

  if (is.null(sweep_result$point_estimates) || is.null(sweep_result$mode)) {
    stop("Input does not appear to be output from sensitivity_sweep().")
  }

  pe <- sweep_result$point_estimates
  rn <- dimnames(pe)[[1]]
  cn <- dimnames(pe)[[2]]

  # Parse cell specification
  parsed <- parse_cell(cell, rn, cn)
  row_idx <- parsed$row
  col_idx <- parsed$col
  cell_label <- paste0(row_idx, " -> ", col_idx)

  if (sweep_result$mode == "1D") {
    # --- 1D: line plot ---
    vals <- pe[row_idx, col_idx, ]
    lam_labels <- dimnames(pe)[[3]]
    lam_vals <- as.numeric(sub("^lambda_", "", lam_labels))

    df <- data.frame(lambda = lam_vals, estimate = vals)

    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$lambda, y = .data$estimate)) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::scale_x_log10() +
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        labels = scales::percent_format()
      ) +
      ggplot2::labs(
        x = expression(lambda ~ "(log scale)"),
        y = "Estimated Support",
        title = title %||% paste("Sensitivity:", cell_label)
      ) +
      ggplot2::theme_minimal(base_size = 12)
  } else {
    # --- 2D: heatmap ---
    mat <- pe[row_idx, col_idx, , ]
    l1_labels <- dimnames(pe)[[3]]
    l2_labels <- dimnames(pe)[[4]]
    l1_vals <- as.numeric(sub("^l1_", "", l1_labels))
    l2_vals <- as.numeric(sub("^l2_", "", l2_labels))

    df <- expand.grid(
      lambda1 = l1_vals, lambda2 = l2_vals,
      KEEP.OUT.ATTRS = FALSE
    )
    df$estimate <- as.vector(mat)
    df$lambda1_f <- factor(df$lambda1)
    df$lambda2_f <- factor(df$lambda2)

    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = .data$lambda2_f,
      y = .data$lambda1_f,
      fill = .data$estimate
    )) +
      ggplot2::geom_tile(color = "white", linewidth = 0.5) +
      ggplot2::scale_fill_viridis_c(
        option = "viridis", name = "Estimate",
        labels = scales::percent_format()
      )

    if (nrow(df) <= 64) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = sprintf("%.2f", .data$estimate)),
        color = "white", size = 3
      )
    }

    p <- p +
      ggplot2::labs(
        x = expression(lambda[2]),
        y = expression(lambda[1]),
        title = title %||% paste("Sensitivity:", cell_label)
      ) +
      ggplot2::theme_minimal(base_size = 12) +
      ggplot2::theme(panel.grid = ggplot2::element_blank())
  }

  p
}


#' Plot Method Comparison
#'
#' Create a dot plot comparing ecological inference estimates across methods,
#' using output from [compare_methods()].
#'
#' @param comparison_result Output of [compare_methods()], a data frame with
#'   a `method` column and one column per group-candidate cell.
#' @param title Optional plot title.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' truth <- matrix(c(0.80, 0.20, 0.15, 0.85),
#'   nrow = 2, byrow = TRUE,
#'   dimnames = list(
#'     c("white", "black"),
#'     c("cand_A", "cand_B")
#'   )
#' )
#' dat <- simulate_election(n_precincts = 40, true_support = truth, seed = 1)
#' comp <- compare_methods(dat,
#'   calibrated_lambda1 = 0.5,
#'   calibrated_lambda2 = 0.25, truth = truth,
#'   sample = 10000, burnin = 2000, thin = 5
#' )
#' plot_comparison(comp)
#' }
#'
#' @export
plot_comparison <- function(comparison_result, title = NULL) {
  check_ggplot2()

  if (!"method" %in% names(comparison_result)) {
    stop(
      "Input must be a data frame with a 'method' column ",
      "(output from compare_methods())."
    )
  }

  # Reshape to long format using base R
  cell_cols <- setdiff(names(comparison_result), "method")
  long <- data.frame(
    method = rep(comparison_result$method, each = length(cell_cols)),
    cell = rep(cell_cols, times = nrow(comparison_result)),
    estimate = as.vector(t(as.matrix(comparison_result[, cell_cols, drop = FALSE]))),
    stringsAsFactors = FALSE
  )

  # Make cell labels more readable
  long$cell <- gsub("_", " -> ", long$cell)

  # Extract truth if present for reference lines
  truth_rows <- long[long$method == "Truth", ]
  est_rows <- long[long$method != "Truth", ]

  # Order methods: put Truth last
  meth_levels <- unique(comparison_result$method)
  meth_levels <- c(
    setdiff(meth_levels, "Truth"),
    intersect(meth_levels, "Truth")
  )
  long$method <- factor(long$method, levels = rev(meth_levels))

  p <- ggplot2::ggplot(long, ggplot2::aes(
    x = .data$estimate,
    y = .data$method
  )) +
    ggplot2::geom_point(size = 3, ggplot2::aes(color = .data$method)) +
    ggplot2::facet_wrap(~cell, scales = "free_x") +
    ggplot2::scale_x_continuous(labels = scales::percent_format()) +
    ggplot2::labs(
      x = "Estimated Support",
      y = NULL,
      title = title %||% "Method Comparison: Group-Conditional Support",
      color = "Method"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none")

  # Add truth reference lines if available
  if (nrow(truth_rows) > 0) {
    p <- p + ggplot2::geom_vline(
      data = truth_rows,
      ggplot2::aes(xintercept = .data$estimate),
      linetype = "dashed", color = "gray40", alpha = 0.7
    )
  }

  p
}


#' Plot Leave-One-Out Cross-Validation Results
#'
#' Visualize LOO-CV results from [loo_cv()]. Two display types are available:
#' a summary bar chart of RMSE and MAE by method, and a per-precinct error
#' dot plot.
#'
#' @param loocv_result Output of [loo_cv()].
#' @param type One of `"summary"` (bar chart of RMSE/MAE) or `"precinct"`
#'   (dot plot of per-precinct errors). Default `"summary"`.
#' @param title Optional plot title.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' # --- Fast Goodman-only example ---
#' dat <- simulate_election(n_precincts = 30, seed = 42)
#' cv <- loo_cv(dat, methods = "goodman", verbose = FALSE)
#' plot_loocv(cv, type = "summary")
#' plot_loocv(cv, type = "precinct")
#'
#' \dontrun{
#' # --- Compare methods ---
#' cv2 <- loo_cv(dat,
#'   methods = c("goodman", "rxc_default"),
#'   sample = 10000, burnin = 2000, thin = 5
#' )
#' plot_loocv(cv2, type = "summary")
#' plot_loocv(cv2, type = "precinct")
#' }
#'
#' @export
plot_loocv <- function(loocv_result, type = c("summary", "precinct"),
                       title = NULL) {
  check_ggplot2()
  type <- match.arg(type)

  if (is.null(loocv_result$summary) || is.null(loocv_result$precinct_errors)) {
    stop("Input does not appear to be output from loo_cv().")
  }

  if (type == "summary") {
    # --- Bar chart of RMSE and MAE ---
    sm <- loocv_result$summary

    # Reshape to long
    long <- rbind(
      data.frame(
        method = sm$method, metric = "RMSE", value = sm$rmse,
        stringsAsFactors = FALSE
      ),
      data.frame(
        method = sm$method, metric = "MAE", value = sm$mae,
        stringsAsFactors = FALSE
      )
    )

    p <- ggplot2::ggplot(long, ggplot2::aes(
      x = .data$method,
      y = .data$value,
      fill = .data$metric
    )) +
      ggplot2::geom_col(position = "dodge", width = 0.7) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.4f", .data$value)),
        position = ggplot2::position_dodge(width = 0.7),
        vjust = -0.3, size = 3
      ) +
      ggplot2::scale_fill_manual(values = c(
        "RMSE" = "steelblue",
        "MAE" = "coral"
      )) +
      ggplot2::labs(
        x = NULL,
        y = "Error",
        fill = "Metric",
        title = title %||% "LOO-CV: Prediction Accuracy by Method"
      ) +
      ggplot2::theme_minimal(base_size = 12)
  } else {
    # --- Per-precinct error dot plot ---
    pe <- loocv_result$precinct_errors
    pe <- pe[!is.na(pe$error), ]

    p <- ggplot2::ggplot(pe, ggplot2::aes(
      x = .data$precinct,
      y = .data$error,
      color = .data$method
    )) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
      ggplot2::geom_point(
        size = 2.5, alpha = 0.8,
        position = ggplot2::position_dodge(width = 0.4)
      ) +
      ggplot2::labs(
        x = "Precinct",
        y = "Prediction Error (predicted - actual)",
        color = "Method",
        title = title %||% "LOO-CV: Per-Precinct Prediction Errors"
      ) +
      ggplot2::theme_minimal(base_size = 12)
  }

  p
}


# --------------------------------------------------------------------------
#  Internal helpers
# --------------------------------------------------------------------------

#' @noRd
check_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop(
      "Package 'ggplot2' is required for plotting. ",
      "Install it with: install.packages('ggplot2')"
    )
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop(
      "Package 'scales' is required for plotting. ",
      "Install it with: install.packages('scales')"
    )
  }
}

#' @noRd
parse_cell <- function(cell, row_names, col_names) {
  if (length(cell) == 2) {
    row_idx <- cell[1]
    col_idx <- cell[2]
  } else if (length(cell) == 1) {
    # Try matching "group_candidate" against all combinations
    found <- FALSE
    for (r in row_names) {
      for (cc in col_names) {
        if (cell == paste0(r, "_", cc)) {
          row_idx <- r
          col_idx <- cc
          found <- TRUE
          break
        }
      }
      if (found) break
    }
    if (!found) {
      # Build available options for error message
      available <- outer(row_names, col_names, paste, sep = "_")
      stop(
        "Cell '", cell, "' not found. Available cells:\n  ",
        paste(as.vector(available), collapse = ", ")
      )
    }
  } else {
    stop(
      "cell must be a length-1 string ('group_candidate') or ",
      "length-2 vector c('group', 'candidate')."
    )
  }

  if (!row_idx %in% row_names) {
    stop(
      "Row '", row_idx, "' not found. Available: ",
      paste(row_names, collapse = ", ")
    )
  }
  if (!col_idx %in% col_names) {
    stop(
      "Column '", col_idx, "' not found. Available: ",
      paste(col_names, collapse = ", ")
    )
  }

  list(row = row_idx, col = col_idx)
}

#' Null-coalescing operator
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x
