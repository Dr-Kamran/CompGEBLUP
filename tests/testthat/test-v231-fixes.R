# Tests for v2.3.1 bug fixes
# Covers: matrix operations, multi-effect VC accuracy, ridge consistency

# ===========================================================================
# Matrix Operations
# ===========================================================================

test_that("make_positive_definite with ridge= adds diagonal ridge", {
  A <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  K <- make_positive_definite(A, ridge = 0.01)

  # Diagonal should increase by exactly 0.01
  expect_equal(diag(K), diag(A) + 0.01, tolerance = 1e-10)

  # Off-diagonal should be unchanged

  expect_equal(K[1, 2], A[1, 2], tolerance = 1e-10)
})

test_that("make_positive_definite with tol= floors eigenvalues (different from ridge)", {
  set.seed(100)
  n <- 20
  X <- matrix(rnorm(n * 10), n, 10)
  K <- X %*% t(X) / 10
  eig <- eigen(K, symmetric = TRUE)
  eig$values[n] <- -0.01
  K_bad <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
  K_bad <- (K_bad + t(K_bad)) / 2

  # tol= floors eigenvalues — changes off-diagonals
  K_floored <- make_positive_definite(K_bad, tol = 0.05, use.cpp = FALSE)
  expect_false(isTRUE(all.equal(K_bad, K_floored, tolerance = 1e-6)))

  # ridge= adds diagonal — preserves off-diagonals
  K_ridged <- make_positive_definite(K_bad, ridge = 0.05)
  expect_equal(K_ridged[1, 2], K_bad[1, 2], tolerance = 1e-10)
})

test_that("add_ridge matches make_positive_definite(ridge=) for PD input", {
  A <- matrix(c(1, 0.5, 0.3, 0.5, 1, 0.4, 0.3, 0.4, 1), 3, 3)
  expect_equal(add_ridge(A, 0.01), make_positive_definite(A, ridge = 0.01),
               tolerance = 1e-10)
})

test_that("build_E_matrix A#A matches Hadamard of A", {
  set.seed(42)
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 100)
  K <- calc_relationship_matrices(gdata, matrices = "A")
  A <- K$A

  AA_built <- build_E_matrix(M = gdata@genotypes, type = "A#A", A = A,
                             min.MAF = 0, use.cpp = TRUE)
  expect_equal(as.matrix(AA_built), as.matrix(A * A), tolerance = 1e-8)
})

test_that("build_GE_matrix produces correct Kronecker dimensions", {
  n <- 20; n_env <- 3
  A <- matrix(runif(n * n), n, n)
  A <- (A + t(A)) / 2; diag(A) <- diag(A) + 1
  rownames(A) <- colnames(A) <- paste0("ID", 1:n)

  K_AE <- build_GE_matrix(A, n_env = n_env)
  expect_equal(nrow(K_AE), n * n_env)
  expect_equal(ncol(K_AE), n * n_env)
})


# ===========================================================================
# Multi-Effect VC Accuracy
# ===========================================================================

test_that("Henderson EM-REML and AI-REML agree for multi-effect model", {
  set.seed(12345)
  gdata <- simulate_genotypes(n_ind = 40, n_snp = 100)
  gdata <- simulate_phenotypes(gdata, n_env = 3, h2 = 0.5, n_qtl = 10)

  K <- calc_relationship_matrices(gdata, matrices = "A")

  design <- build_design_matrices(
    pheno = gdata@phenotypes, K_matrices = K, effects = c("A", "ENV")
  )

  y <- as.matrix(gdata@phenotypes$trait)
  X <- design$X
  Z_list <- design$Z_list
  K_list <- lapply(design$K_list, function(K) K + diag(0.001, nrow(K)))

  r_h <- CompGEBLUP:::fit_mme_henderson(
    y = y, X = X, Z_list = Z_list, K_list = K_list,
    nIters = 500, tolParConvLL = 1e-6, tolParInv = 1e-6, verbose = FALSE)

  r_ai <- CompGEBLUP:::fit_mme_henderson_ai(
    y = y, X = X, Z_list = Z_list, K_list = K_list,
    nIters = 500, tolParConvLL = 1e-6, tolParInv = 1e-6, verbose = FALSE)

  expect_true(r_h$convergence)
  expect_true(r_ai$convergence)

  # VCs should agree closely
  for (j in seq_along(r_h$theta)) {
    expect_equal(unname(r_h$theta[j]), unname(r_ai$theta[j]),
                 tolerance = 0.15, label = paste("VC", j))
  }

  # BLUPs highly correlated
  for (i in seq_along(r_h$u)) {
    expect_true(cor(as.vector(r_h$u[[i]]), as.vector(r_ai$u[[i]])) > 0.95)
  }
})

test_that("Full EM-REML sigma2_e produces reasonable VCs for multi-effect", {
  set.seed(999)
  gdata <- simulate_genotypes(n_ind = 50, n_snp = 200)
  gdata <- simulate_phenotypes(gdata, n_env = 4, h2 = 0.3, n_qtl = 10)

  K <- calc_relationship_matrices(gdata, matrices = "A")

  model <- fit_gblup(gdata, K, effects = c("A", "ENV"),
                     ridge = 0.001, nIters = 200,
                     verbose = FALSE, use.emma = FALSE)

  vc <- get_varcomp(model)
  expect_true(all(vc$Variance > 0))
  h2 <- vc$Variance[vc$Component == "A"] / sum(vc$Variance)
  expect_true(h2 > 0.01 && h2 < 0.95)
})


# ===========================================================================
# Ridge Consistency
# ===========================================================================

test_that("ridge_per_matrix in fit_gblup applies per-component ridge", {
  set.seed(123)
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 5)

  K <- calc_relationship_matrices(gdata, matrices = c("A", "D"))

  model <- fit_gblup(gdata, K, effects = c("A", "D"),
                     ridge = 0.001,
                     ridge_per_matrix = c(A = 0.001, D = 0.05),
                     nIters = 200, verbose = FALSE, use.emma = FALSE)

  expect_s4_class(model, "GBLUPModel")
  expect_true(all(get_varcomp(model)$Variance >= 0))
})

test_that("cv_gblup with recompute_K=TRUE uses diagonal ridge", {
  set.seed(42)
  gdata <- simulate_genotypes(n_ind = 40, n_snp = 100)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)

  K <- calc_relationship_matrices(gdata, matrices = "A")

  cv_res <- cv_gblup(gdata, K, effects = "A",
                     scheme = "CV1", n_folds = 3, n_reps = 1,
                     ridge = 0.001,
                     ridge_per_matrix = c(A = 0.001),
                     recompute_K = TRUE, verbose = FALSE)

  expect_s4_class(cv_res, "CVResult")
  expect_true(all(!is.na(cv_res@metrics$predictive_ability)))
})

test_that("tolParInv decoupled from user ridge in fit_gblup", {
  set.seed(42)
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.5, n_qtl = 5)

  K <- calc_relationship_matrices(gdata, matrices = "A")

  # Should not error with large ridge (tolParInv now fixed at 1e-6)
  model <- fit_gblup(gdata, K, effects = "A",
                     ridge = 0.05, nIters = 100,
                     verbose = FALSE, use.emma = FALSE)
  expect_s4_class(model, "GBLUPModel")
})
