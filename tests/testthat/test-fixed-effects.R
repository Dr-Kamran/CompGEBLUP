# Test Fixed Effects Functionality

test_that("fit_gblup accepts fixed parameter with formula", {
  set.seed(123)
  
  # Generate data
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  # Add TESTER column to phenotypes
  gdata@phenotypes$TESTER <- rep(c("T1", "T2"), length.out = nrow(gdata@phenotypes))
  
  # Calculate matrices
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Fit model with fixed effects
  expect_silent(
    model <- fit_gblup(gdata, K, effects = "A", fixed = ~ TESTER, verbose = FALSE)
  )
  
  # Tests
  expect_s4_class(model, "GBLUPModel")
  expect_true(!is.null(model@model))
  expect_true(is.data.frame(model@gebv))
  expect_true(is.data.frame(model@varcomp))
})

test_that("fit_gblup with fixed effects produces reasonable results", {
  set.seed(123)
  
  # Generate data
  gdata <- simulate_genotypes(n_ind = 40, n_snp = 100)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  # Add TESTER column with strong effect
  gdata@phenotypes$TESTER <- rep(c("T1", "T2"), length.out = nrow(gdata@phenotypes))
  # Add tester effect to trait
  tester_effect <- ifelse(gdata@phenotypes$TESTER == "T1", 2, -2)
  gdata@phenotypes$trait <- gdata@phenotypes$trait + tester_effect
  
  # Calculate matrices
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Fit model WITHOUT fixed effects
  model_no_fixed <- fit_gblup(gdata, K, effects = "A", verbose = FALSE)
  
  # Fit model WITH fixed effects
  model_with_fixed <- fit_gblup(gdata, K, effects = "A", fixed = ~ TESTER, verbose = FALSE)
  
  # The model with fixed effects should have non-zero genetic variance
  # while the model without fixed effects might have low genetic variance
  # because tester effects are absorbed into residual
  var_a_no_fixed <- model_no_fixed@varcomp$Variance[model_no_fixed@varcomp$Component == "A"]
  var_a_with_fixed <- model_with_fixed@varcomp$Variance[model_with_fixed@varcomp$Component == "A"]
  
  # Both should have valid variance components
  expect_true(var_a_no_fixed >= 0)
  expect_true(var_a_with_fixed >= 0)
  
  # Both should produce GEBVs with meaningful variation (not all zeros)
  gebv_var_no_fixed <- var(model_no_fixed@gebv$GEBV, na.rm = TRUE)
  gebv_var_with_fixed <- var(model_with_fixed@gebv$GEBV, na.rm = TRUE)
  expect_true(gebv_var_no_fixed > 1e-6)
  expect_true(gebv_var_with_fixed > 1e-6)
})

test_that("fit_gblup accepts NULL for fixed parameter (default)", {
  set.seed(123)
  
  # Generate data
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  # Calculate matrices
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Fit model with explicit NULL (intercept only)
  expect_silent(
    model <- fit_gblup(gdata, K, effects = "A", fixed = NULL, verbose = FALSE)
  )
  
  # Tests
  expect_s4_class(model, "GBLUPModel")
  expect_true(!is.null(model@model))
})

test_that("fit_gblup accepts pre-built design matrix for fixed", {
  set.seed(123)
  
  # Generate data
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  # Create a design matrix manually
  n <- nrow(gdata@phenotypes)
  X_fixed <- matrix(c(rep(1, n), rnorm(n)), ncol = 2)
  colnames(X_fixed) <- c("(Intercept)", "covariate")
  
  # Calculate matrices
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Fit model with pre-built matrix
  expect_silent(
    model <- fit_gblup(gdata, K, effects = "A", fixed = X_fixed, verbose = FALSE)
  )
  
  # Tests
  expect_s4_class(model, "GBLUPModel")
  expect_true(!is.null(model@model))
})

test_that("cv_gblup accepts fixed parameter with formula", {
  set.seed(123)
  
  # Generate data
  gdata <- simulate_genotypes(n_ind = 40, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  # Add TESTER column to phenotypes
  gdata@phenotypes$TESTER <- rep(c("T1", "T2"), length.out = nrow(gdata@phenotypes))
  
  # Calculate matrices
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Run CV with fixed effects
  expect_silent(
    cv_res <- cv_gblup(gdata, K, effects = "A", fixed = ~ TESTER, 
                       scheme = "CV1", n_folds = 3, verbose = FALSE)
  )
  
  # Tests
  expect_s4_class(cv_res, "CVResult")
  expect_true(is.data.frame(cv_res@metrics))
  expect_true(is.data.frame(cv_res@predictions))
})

test_that("cv_gblup with fixed effects produces valid predictions", {
  set.seed(123)
  
  # Generate data
  gdata <- simulate_genotypes(n_ind = 50, n_snp = 100)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  # Add TESTER column
  gdata@phenotypes$TESTER <- rep(c("T1", "T2"), length.out = nrow(gdata@phenotypes))
  
  # Calculate matrices
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Run CV with fixed effects
  cv_res <- cv_gblup(gdata, K, effects = "A", fixed = ~ TESTER,
                     scheme = "CV1", n_folds = 3, verbose = FALSE)
  
  # Check that we got valid predictions
  expect_true(nrow(cv_res@predictions) > 0)
  expect_true("predicted" %in% names(cv_res@predictions))
  expect_true("observed" %in% names(cv_res@predictions))
  
  # Check metrics
  expect_true(nrow(cv_res@metrics) == 3)
  expect_true(all(c("predictive_ability", "mse") %in% names(cv_res@metrics)))
  
  # Predictive ability should be a valid correlation
  mean_pa <- mean(cv_res@metrics$predictive_ability, na.rm = TRUE)
  expect_true(!is.na(mean_pa))
  expect_true(abs(mean_pa) <= 1)
})

test_that("cv_gblup accepts NULL for fixed parameter (default)", {
  set.seed(123)
  
  # Generate data
  gdata <- simulate_genotypes(n_ind = 40, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  # Calculate matrices
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Run CV with explicit NULL
  expect_silent(
    cv_res <- cv_gblup(gdata, K, effects = "A", fixed = NULL,
                       scheme = "CV1", n_folds = 3, verbose = FALSE)
  )
  
  # Tests
  expect_s4_class(cv_res, "CVResult")
})

test_that("fixed effects work with multiple covariates", {
  set.seed(123)
  
  # Generate data
  gdata <- simulate_genotypes(n_ind = 40, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  # Add multiple covariates
  gdata@phenotypes$TESTER <- rep(c("T1", "T2", "T3"), length.out = nrow(gdata@phenotypes))
  gdata@phenotypes$BLOCK <- rep(c("B1", "B2"), length.out = nrow(gdata@phenotypes))
  
  # Calculate matrices
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Fit model with multiple fixed effects
  expect_silent(
    model <- fit_gblup(gdata, K, effects = "A", fixed = ~ TESTER + BLOCK, verbose = FALSE)
  )
  
  # Tests
  expect_s4_class(model, "GBLUPModel")
  expect_true(!is.null(model@model))
})

test_that("cv_fold properly forwards fixed parameter", {
  set.seed(123)
  
  # Generate data
  gdata <- simulate_genotypes(n_ind = 40, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  # Add TESTER column
  gdata@phenotypes$TESTER <- rep(c("T1", "T2"), length.out = nrow(gdata@phenotypes))
  
  # Calculate matrices
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Create a fold manually
  pheno <- gdata@phenotypes
  n <- nrow(pheno)
  test_idx <- 1:10
  train_idx <- 11:n
  
  fold <- list(
    test = test_idx,
    train = train_idx,
    test_individuals = unique(pheno$GID[test_idx])
  )
  
  # Run cv_fold with fixed effects
  # This is an internal function, so we need to use ::: to access it
  expect_silent(
    result <- CompGEBLUP:::cv_fold(gdata, K, effects = "A", fold = fold, 
                                    ridge = 0.001, recompute_K = FALSE, 
                                    verbose = FALSE, fixed = ~ TESTER)
  )
  
  # Check that result is valid
  expect_true(is.list(result))
  expect_true("observed" %in% names(result))
  expect_true("predicted" %in% names(result))
})
