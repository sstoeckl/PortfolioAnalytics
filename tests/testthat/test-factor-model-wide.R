###############################################################################
# tests/testthat/test-factor-model-wide.R
#
# statistical.factor.model() is a principal component decomposition, and PCA
# does not require more observations than assets. The function nonetheless
# refused whenever m < N, which is the one case a factor structure is normally
# reached for: a universe wider than its history, where the sample covariance
# is singular.
#
# What the fit actually needs is k usable components. prcomp() returns
# min(m, N) of them and min(m - 1, N) carry variance once the data is centred,
# so that is the bound these tests pin.
###############################################################################

skip_on_cran()
library(PortfolioAnalytics)

utils::data(edhec)
R <- edhec[, 1:8]                       # 8 assets, plenty of months
wide <- R[1:5, ]                        # 5 observations, 8 assets: m < N

test_that("a fit with fewer observations than assets is accepted", {
  expect_no_error(statistical.factor.model(wide, k = 2))
})

test_that("the wide fit produces a covariance of the right shape and rank", {
  S <- extractCovariance(statistical.factor.model(wide, k = 2))
  expect_equal(dim(S), c(ncol(wide), ncol(wide)))
  expect_equal(qr(S)$rank, ncol(wide))
  expect_gt(min(eigen(S, symmetric = TRUE, only.values = TRUE)$values), 0)
})

test_that("k is bounded by the components the data supports, not by N", {
  m <- nrow(wide)
  # min(m - 2, N) = 3 here. extractCovariance() divides by m - k - 1, so
  # k = m - 1 would return an infinite covariance instead of an error.
  expect_no_error(statistical.factor.model(wide, k = m - 2L))
  expect_error(statistical.factor.model(wide, k = m - 1L),
               "more factors than the data supports")
  expect_error(statistical.factor.model(wide, k = 99L),
               "more factors than the data supports")
})

test_that("the asset count is the other half of the bound", {
  # The bound is min(m - 2, N). The wide fixture above never reaches the N
  # half of it: with m = 5 the observation term binds at every k, so a
  # regression that dropped N from the bound would go unnoticed there.
  # A tall fixture puts N in charge.
  tall <- R[1:60, ]                     # m = 60, N = 8, so m - 2 = 58 > N
  N    <- ncol(tall)

  expect_no_error(statistical.factor.model(tall, k = N))
  expect_error(statistical.factor.model(tall, k = N + 1L),
               "more factors than the data supports")

  # k = N uses every component, so the residuals vanish and the factor
  # covariance is the sample covariance. That is what makes N the ceiling.
  S <- extractCovariance(statistical.factor.model(tall, k = N))
  expect_equal(as.numeric(S),
               as.numeric(stats::cov(zoo::coredata(tall))),
               tolerance = 1e-12)
})

test_that("the square case k = N = m is refused, and that is a change", {
  # master ran this one: a full basis reconstructs the data, the residuals are
  # zero, and the negative denominator m - k - 1 = -1 has nothing to divide,
  # so extractCovariance() came back with the sample covariance of a square
  # panel -- finite, and singular, which is the case a factor model exists to
  # avoid. It is refused here deliberately. Admitting it would make the
  # admissible set stop being an interval: one more observation, m = N + 1
  # with k = N, divides by zero instead.
  square <- R[seq_len(ncol(R)), ]       # m = N = 8
  expect_error(statistical.factor.model(square, k = ncol(square)),
               "requests more factors than the data supports")
  expect_error(statistical.factor.model(R[seq_len(ncol(R) + 1L), ], k = ncol(R)),
               "requests more factors than the data supports")
  expect_no_error(statistical.factor.model(square, k = ncol(square) - 2L))
})

test_that("every accepted k yields a finite covariance", {
  # The bound exists for this reason, so check the reason rather than
  # only the boundary.
  for (k in seq_len(nrow(wide) - 2L)) {
    S <- extractCovariance(statistical.factor.model(wide, k = k))
    expect_true(all(is.finite(S)), info = paste("k =", k))
    expect_gt(min(diag(S)), 0)
  }
})

test_that("the message says what was asked for and what was available", {
  msg <- tryCatch(statistical.factor.model(wide, k = 99L),
                  error = function(e) conditionMessage(e))
  expect_match(msg, "k = 99")
  expect_match(msg, as.character(nrow(wide)))
  expect_match(msg, as.character(ncol(wide)))
})

test_that("too few observations for any factor are refused", {
  # m - 2 < 1 leaves no admissible k at all.
  expect_error(statistical.factor.model(R[1, ], k = 1),
               "at least three observations")
  expect_error(statistical.factor.model(R[1:2, ], k = 1),
               "at least three observations")
  expect_no_error(statistical.factor.model(R[1:3, ], k = 1))
})

test_that("k <= 0 is still refused, and before the new check", {
  expect_error(statistical.factor.model(wide, k = 0), "positive integer")
  expect_error(statistical.factor.model(wide, k = -1), "positive integer")
})

test_that("the tall case is unchanged", {
  # Where the old guard allowed a fit, the fit must be the same one.
  tall <- R[1:60, ]
  fit  <- statistical.factor.model(tall, k = 3)
  S    <- extractCovariance(fit)

  # Reproduce the decomposition independently.
  x     <- zoo::coredata(tall)
  pc    <- prcomp(x)
  betas <- pc$rotation[, 1:3, drop = FALSE]
  f     <- x %*% betas
  resid <- x - f %*% t(betas)
  resid <- sweep(resid, 2, colMeans(resid))
  # extractCovariance() divides the residual sums of squares by m - k - 1,
  # not by m - 1, so var() is not the right denominator here.
  stock_m2 <- colSums(resid^2) / (nrow(x) - 3L - 1L)
  S_ref <- betas %*% stats::cov(f) %*% t(betas) + diag(stock_m2)

  expect_equal(as.numeric(S), as.numeric(S_ref), tolerance = 1e-8)
})
