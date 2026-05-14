#' Simulate a Two-Group, Two-Candidate Election with Known Truth
#'
#' Generates a precinct-level dataset where the true group-conditional vote
#' shares are known. Useful for benchmarking ecological inference methods and
#' for the calibration workflow in [calibrate_rxc()].
#'
#' @param n_precincts Integer. Number of precincts to simulate.
#' @param true_support A 2x2 matrix of true group-conditional vote shares, with
#'   rows summing to 1. Rows are demographic groups (default: white, black),
#'   columns are candidates (default: cand_A, cand_B).
#' @param precinct_size_mean Mean precinct size (voters per precinct).
#' @param precinct_size_sd Standard deviation of precinct sizes.
#' @param noise Standard deviation of precinct-level noise added to the true
#'   group rates before drawing votes. Larger values produce noisier data.
#' @param composition_shape Two-element vector giving the shape parameters for
#'   the Beta distribution generating precinct racial composition. Default
#'   `c(2, 2)` gives mixed precincts; `c(0.5, 0.5)` gives more polarized
#'   (U-shaped) composition; `c(5, 5)` gives more uniform composition.
#' @param seed Optional random seed for reproducibility.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{precinct}{Precinct identifier.}
#'   \item{total}{Total voters in the precinct.}
#'   \item{n_white, n_black}{Voter counts by group.}
#'   \item{pct_white, pct_black}{Voter proportions by group.}
#'   \item{cand_A, cand_B}{Vote counts by candidate.}
#'   \item{pct_A, pct_B}{Vote proportions by candidate.}
#' }
#'
#' @examples
#' # --- Basic usage with default truth (80-20 / 15-85) ---
#' dat <- simulate_election(n_precincts = 40, seed = 1)
#' head(dat)
#'
#' # Inspect the embedded ground truth
#' attr(dat, "true_support")
#'
#' # --- Custom truth matrix ---
#' truth <- matrix(c(0.70, 0.30, 0.25, 0.75), nrow = 2, byrow = TRUE,
#'                 dimnames = list(c("white", "black"),
#'                                 c("cand_A", "cand_B")))
#' dat2 <- simulate_election(n_precincts = 50, true_support = truth, seed = 2)
#'
#' # --- Polarized racial composition (U-shaped: mostly segregated precincts) ---
#' dat_polar <- simulate_election(n_precincts = 40,
#'                                composition_shape = c(0.5, 0.5),
#'                                seed = 3)
#' summary(dat_polar$pct_white)  # most precincts near 0 or 1
#'
#' # --- Homogeneous composition (all precincts ~50/50) ---
#' # This creates a hard identification problem for EI methods
#' dat_homo <- simulate_election(n_precincts = 40,
#'                               composition_shape = c(50, 50),
#'                               seed = 4)
#' summary(dat_homo$pct_white)  # all precincts near 0.50
#'
#' # --- Higher noise (noisier precinct-level rates) ---
#' dat_noisy <- simulate_election(n_precincts = 40, noise = 0.15, seed = 5)
#'
#' @export
simulate_election <- function(n_precincts = 50,
                              true_support = matrix(
                                c(0.80, 0.20,
                                  0.15, 0.85),
                                nrow = 2, byrow = TRUE,
                                dimnames = list(c("white", "black"),
                                                c("cand_A", "cand_B"))),
                              precinct_size_mean = 800,
                              precinct_size_sd = 100,
                              noise = 0.05,
                              composition_shape = c(2, 2),
                              seed = NULL) {

  if (!is.null(seed)) set.seed(seed)

  stopifnot(is.matrix(true_support),
            nrow(true_support) == 2,
            ncol(true_support) == 2)
  if (any(abs(rowSums(true_support) - 1) > 1e-8)) {
    warning("Rows of true_support do not sum to 1; rescaling.")
    true_support <- true_support / rowSums(true_support)
  }

  row_names <- rownames(true_support) %||% c("white", "black")
  col_names <- colnames(true_support) %||% c("cand_A", "cand_B")

  pct_white <- stats::rbeta(n_precincts, composition_shape[1], composition_shape[2])
  pct_black <- 1 - pct_white
  total <- round(stats::rnorm(n_precincts, precinct_size_mean, precinct_size_sd))
  total[total < 100] <- 100

  n_white <- round(pct_white * total)
  n_black <- total - n_white

  cand_A <- integer(n_precincts)
  cand_B <- integer(n_precincts)

  for (i in seq_len(n_precincts)) {
    w_a_rate <- clamp(true_support[1, 1] + stats::rnorm(1, 0, noise), 0.01, 0.99)
    b_a_rate <- clamp(true_support[2, 1] + stats::rnorm(1, 0, noise), 0.01, 0.99)

    w_a <- stats::rbinom(1, n_white[i], w_a_rate)
    b_a <- stats::rbinom(1, n_black[i], b_a_rate)

    cand_A[i] <- w_a + b_a
    cand_B[i] <- total[i] - cand_A[i]
  }

  out <- data.frame(
    precinct  = seq_len(n_precincts),
    total     = total,
    n_white   = n_white,
    n_black   = n_black,
    pct_white = n_white / total,
    pct_black = n_black / total,
    cand_A    = cand_A,
    cand_B    = cand_B,
    pct_A     = cand_A / total,
    pct_B     = cand_B / total
  )

  # Stash truth and naming for downstream use
  attr(out, "true_support") <- true_support
  attr(out, "row_names")    <- row_names
  attr(out, "col_names")    <- col_names

  out
}

# Internal helper: clamp x to [lo, hi]
clamp <- function(x, lo, hi) pmin(pmax(x, lo), hi)

# Internal helper: null-coalescing operator
`%||%` <- function(a, b) if (is.null(a)) b else a
