test_that("simulate_election returns a properly shaped data frame", {
  dat <- simulate_election(n_precincts = 20, seed = 1)
  expect_s3_class(dat, "data.frame")
  expect_equal(nrow(dat), 20)
  expect_true(all(c("total", "n_white", "n_black",
                    "cand_A", "cand_B", "pct_white", "pct_black",
                    "pct_A", "pct_B") %in% names(dat)))
})

test_that("simulate_election row counts add up", {
  dat <- simulate_election(n_precincts = 50, seed = 2)
  expect_equal(dat$n_white + dat$n_black, dat$total)
  expect_equal(dat$cand_A + dat$cand_B, dat$total)
})

test_that("simulate_election respects true_support direction", {
  truth <- matrix(c(0.90, 0.10, 0.10, 0.90), nrow = 2, byrow = TRUE,
                  dimnames = list(c("white", "black"),
                                  c("cand_A", "cand_B")))
  dat <- simulate_election(n_precincts = 200, true_support = truth,
                           noise = 0.01, seed = 3)
  # In highly-white precincts, cand_A should dominate
  white_heavy <- dat[dat$pct_white > 0.8, ]
  black_heavy <- dat[dat$pct_black > 0.8, ]
  expect_gt(mean(white_heavy$pct_A), 0.7)
  expect_gt(mean(black_heavy$pct_B), 0.7)
})

test_that("simulate_election normalizes truth rows that don't sum to 1", {
  bad_truth <- matrix(c(0.8, 0.4, 0.2, 0.6), nrow = 2, byrow = TRUE,
                      dimnames = list(c("white", "black"),
                                      c("cand_A", "cand_B")))
  expect_warning(simulate_election(n_precincts = 10, true_support = bad_truth,
                                   seed = 4),
                 regexp = "rescal")
})
