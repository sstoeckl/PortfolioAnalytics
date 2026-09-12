###############################################################################
# tests/testthat/test-stat-factor-model.R
#
# Tests for R/stat.factor.model.R
#
# Covers:
#   statistical.factor.model()  — PCA-based statistical factor model
#   center()                    — column-centering of a matrix
#   extractCovariance()         — covariance matrix from stfm object
#   extractCoskewness()         — coskewness matrix from stfm object
#   extractCokurtosis()         — cokurtosis matrix from stfm object
#
# Copyright (c) 2004-2026 Brian G. Peterson, Peter Carl, Ross Bennett
# License: GPL-3
###############################################################################

library(PortfolioAnalytics)
library(PerformanceAnalytics)

skip_on_cran()

# ---------------------------------------------------------------------------
# File-scope setup — runs once before any test_that() block.
# edhec5 (293 obs × 5 assets) and make_mines_portf() are provided by
# helper-portfolioanalytics.R and do NOT need to be redefined here.
# ---------------------------------------------------------------------------

N <- ncol(edhec5)   # 5 assets
m <- nrow(edhec5)   # 293 observations

# Single-factor and two-factor models fitted once for all tests
sfm_k1 <- statistical.factor.model(edhec5, k = 1)
sfm_k2 <- statistical.factor.model(edhec5, k = 2)

# Plain numeric matrix from edhec5 for center() tests
x_mat <- coredata(edhec5)


# ===========================================================================
# statistical.factor.model() — structure tests
# ===========================================================================

test_that("statistical.factor.model(k=1) returns an object of class 'stfm'", {
  expect_s3_class(sfm_k1, "stfm")
})

test_that("stfm object contains all required named slots", {
  expect_true(!is.null(sfm_k1$factor_loadings),    label = "slot factor_loadings exists")
  expect_true(!is.null(sfm_k1$factor_realizations), label = "slot factor_realizations exists")
  expect_true(!is.null(sfm_k1$residuals),           label = "slot residuals exists")
  expect_true(!is.null(sfm_k1$m),                   label = "slot m exists")
  expect_true(!is.null(sfm_k1$k),                   label = "slot k exists")
  expect_true(!is.null(sfm_k1$N),                   label = "slot N exists")
})

test_that("stfm scalar slots m, k, N record correct values for k=1", {
  expect_equal(sfm_k1$m, m)
  expect_equal(sfm_k1$k, 1L)
  expect_equal(sfm_k1$N, N)
})

test_that("factor_loadings for k=1 is a numeric vector of length N (R drops dim on single-column index)", {
  fl <- sfm_k1$factor_loadings
  # prcomp rotation[, 1] returns a named numeric vector, not a matrix
  expect_true(is.numeric(fl),     label = "factor_loadings is numeric")
  expect_false(is.matrix(fl),     label = "factor_loadings is a vector (not a matrix) when k=1")
  expect_equal(length(fl), N,     label = "factor_loadings has length N")
})

test_that("factor_realizations is an m x 1 matrix for k=1", {
  fr <- sfm_k1$factor_realizations
  expect_true(is.matrix(fr),      label = "factor_realizations is a matrix")
  expect_equal(nrow(fr), m,       label = "factor_realizations has m rows")
  expect_equal(ncol(fr), 1L,      label = "factor_realizations has 1 column")
})

test_that("residuals is an m x N matrix for k=1", {
  res <- sfm_k1$residuals
  expect_true(is.matrix(res),     label = "residuals is a matrix")
  expect_equal(nrow(res), m,      label = "residuals has m rows")
  expect_equal(ncol(res), N,      label = "residuals has N columns")
})

test_that("statistical.factor.model(k=2) produces a two-factor model with correct dimensions", {
  expect_equal(sfm_k2$k, 2L,
    label = "k slot equals 2L")
  expect_equal(ncol(sfm_k2$factor_loadings), 2L,
    label = "factor_loadings has 2 columns")
  expect_equal(ncol(sfm_k2$factor_realizations), 2L,
    label = "factor_realizations has 2 columns")
  expect_equal(nrow(sfm_k2$factor_loadings), N,
    label = "factor_loadings still has N rows")
  expect_equal(nrow(sfm_k2$residuals), m,
    label = "residuals still has m rows")
})

test_that("statistical.factor.model() accepts fewer observations than assets", {
  # This test asserted the m < N guard that this change removes, so it is
  # rewritten rather than added to: the guard is the whole subject of both.
  #
  # The fit is a principal component decomposition, which needs k usable
  # components and not m >= N. The bound that does apply is
  # k <= min(m - 2, N), because extractCovariance() divides the residual
  # sums of squares by m - k - 1.
  too_few <- edhec5[1:4, ]        # 4 observations, 5 assets: m < N
  expect_no_error(statistical.factor.model(too_few, k = 1))
  expect_error(
    statistical.factor.model(too_few, k = 3),
    "more factors than the data supports"
  )
})


# ===========================================================================
# center() — column-centering function
# ===========================================================================

test_that("center() returns a matrix with all column means equal to zero (tol 1e-10)", {
  centered <- center(x_mat)
  col_means <- colMeans(centered)
  expect_true(
    all(abs(col_means) < 1e-10),
    label = "every column mean is within 1e-10 of zero"
  )
})

test_that("center() preserves the matrix dimensions of its input", {
  centered <- center(x_mat)
  expect_equal(dim(centered), dim(x_mat))
})

test_that("center() errors with an informative message on non-matrix input", {
  # data.frame
  expect_error(center(as.data.frame(x_mat)), "x must be a matrix")
  # plain numeric vector
  expect_error(center(as.numeric(x_mat[, 1])), "x must be a matrix")
})


# ===========================================================================
# extractCovariance() — covariance matrix from stfm
# ===========================================================================

test_that("extractCovariance() returns a symmetric N x N matrix", {
  cov_mat <- extractCovariance(sfm_k1)
  expect_true(is.matrix(cov_mat),
    label = "result is a matrix")
  expect_equal(dim(cov_mat), c(N, N),
    label = "dimensions are N x N")
  expect_true(
    isSymmetric(cov_mat, tol = 1e-10),
    label = "matrix is symmetric"
  )
})

test_that("extractCovariance() result is positive semi-definite (all eigenvalues >= -1e-8)", {
  cov_mat <- extractCovariance(sfm_k1)
  eigs    <- eigen(cov_mat, symmetric = TRUE, only.values = TRUE)$values
  expect_true(
    min(eigs) >= -1e-8,
    label = "minimum eigenvalue is >= -1e-8"
  )
})

test_that("extractCovariance() also works for a two-factor model", {
  cov_mat2 <- extractCovariance(sfm_k2)
  expect_true(is.matrix(cov_mat2))
  expect_equal(dim(cov_mat2), c(N, N))
  expect_true(isSymmetric(cov_mat2, tol = 1e-10))
})


# ===========================================================================
# extractCoskewness() — coskewness matrix from stfm
# ===========================================================================

test_that("extractCoskewness() returns a matrix with N rows and N^2 columns", {
  m3 <- extractCoskewness(sfm_k1)
  expect_true(is.matrix(m3),
    label = "coskewness result is a matrix")
  expect_equal(nrow(m3), N,
    label = "coskewness has N rows")
  expect_equal(ncol(m3), N^2,
    label = "coskewness has N^2 columns")
})

test_that("extractCoskewness() also produces the correct shape for k=2", {
  m3_k2 <- extractCoskewness(sfm_k2)
  expect_equal(dim(m3_k2), c(N, N^2))
})


# ===========================================================================
# extractCokurtosis() — cokurtosis matrix from stfm
# ===========================================================================

test_that("extractCokurtosis() returns a matrix with N rows and N^3 columns", {
  m4 <- extractCokurtosis(sfm_k1)
  expect_true(is.matrix(m4),
    label = "cokurtosis result is a matrix")
  expect_equal(nrow(m4), N,
    label = "cokurtosis has N rows")
  expect_equal(ncol(m4), N^3,
    label = "cokurtosis has N^3 columns")
})

test_that("extractCokurtosis() also produces the correct shape for k=2", {
  m4_k2 <- extractCokurtosis(sfm_k2)
  expect_equal(dim(m4_k2), c(N, N^3))
})


# ===========================================================================
# Integration smoke test — stfm output feeds set.portfolio.moments(boudt)
# ===========================================================================

test_that("set.portfolio.moments(method='boudt') runs without error for an ES portfolio", {
  portf_es      <- make_mines_portf(edhec5)
  moments_boudt <- set.portfolio.moments(edhec5, portf_es, method = "boudt")
  expect_true(is.list(moments_boudt),
    label = "boudt moments is a list")
  expect_true(!is.null(moments_boudt$mu),
    label = "boudt moments contains mu")
  expect_true(!is.null(moments_boudt$sigma),
    label = "boudt moments contains sigma")
})

test_that("set.portfolio.moments(method='boudt') populates higher moments for an ES portfolio", {
  portf_es      <- make_mines_portf(edhec5)
  moments_boudt <- set.portfolio.moments(edhec5, portf_es, method = "boudt")
  expect_true(!is.null(moments_boudt$m3),
    label = "boudt ES moments contains m3")
  expect_true(!is.null(moments_boudt$m4),
    label = "boudt ES moments contains m4")
})
