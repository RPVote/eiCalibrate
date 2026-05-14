test_that("score_estimates computes correct losses", {
  truth <- matrix(c(0.80, 0.20, 0.15, 0.85), nrow = 2, byrow = TRUE,
                  dimnames = list(c("white", "black"),
                                  c("cand_A", "cand_B")))
  est <- truth  # perfect estimate
  s <- score_estimates(est, truth)
  expect_equal(s$rmse, 0)
  expect_equal(s$mae, 0)
  expect_equal(s$max_err, 0)

  est2 <- truth + matrix(c(0.1, -0.1, -0.1, 0.1), 2, 2)
  s2 <- score_estimates(est2, truth)
  expect_equal(s2$mae, 0.1)
  expect_equal(s2$max_err, 0.1)
})

test_that("score_estimates errors on dimension mismatch", {
  expect_error(score_estimates(matrix(0, 2, 2), matrix(0, 3, 3)),
               "same dimensions")
})

test_that("fit_ei goodman returns a 2x2 matrix in [0,1]", {
  dat <- simulate_election(n_precincts = 50, seed = 10)
  fit <- fit_ei(dat, method = "goodman")
  expect_equal(fit$method, "goodman")
  expect_equal(dim(fit$estimates), c(2, 2))
  expect_true(all(fit$estimates >= 0 & fit$estimates <= 1))
})

test_that("fit_ei goodman recovers truth direction on clean data", {
  truth <- matrix(c(0.90, 0.10, 0.10, 0.90), nrow = 2, byrow = TRUE,
                  dimnames = list(c("white", "black"),
                                  c("cand_A", "cand_B")))
  dat <- simulate_election(n_precincts = 200, true_support = truth,
                           noise = 0.01, seed = 11)
  fit <- fit_ei(dat, method = "goodman")
  # White should support A more than black does
  expect_gt(fit$estimates["white", "cand_A"],
            fit$estimates["black", "cand_A"])
})
