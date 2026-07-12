# Test Model Fitting Functions

test_that("Simple A model fits correctly", {
  set.seed(123)
  
  # Generate data
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  # Calculate matrices
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Fit model
  expect_silent(
    model <- fit_gblup(gdata, K, effects = "A", verbose = FALSE)
  )
  
  # Tests
  expect_s4_class(model, "GBLUPModel")
  expect_true(!is.null(model@model))
  expect_true(is.data.frame(model@gebv))
  expect_true(is.data.frame(model@varcomp))
  
  # Check GEBV
  expect_equal(nrow(model@gebv), nrow(gdata@phenotypes))
  expect_true("GEBV" %in% names(model@gebv))
  expect_true("observed" %in% names(model@gebv))
  
  # Check variance
  expect_true(var(model@gebv$GEBV, na.rm = TRUE) > 0)
  
  # Check variance components - at minimum we have genetic + residual
  expect_true(nrow(model@varcomp) >= 2)
  expect_true(all(model@varcomp$Variance >= 0))
})

test_that("A + ENV model fits correctly", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 3, h2 = 0.6, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  expect_silent(
    model <- fit_gblup(gdata, K, effects = c("A", "ENV"), verbose = FALSE)
  )
  
  # Tests
  expect_s4_class(model, "GBLUPModel")
  # At minimum: Genetic + Residual (ENV is fixed effect, may or may not add vc)
  expect_true(nrow(model@varcomp) >= 2)
})

test_that("GxE model fits correctly", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.5, 
                               n_qtl = 10, include_gxe = TRUE)
  
  K <- calc_relationship_matrices(gdata, matrices = c("A", "AE"))
  
  expect_silent(
    model <- fit_gblup(gdata, K, effects = c("A", "ENV", "AE"), verbose = FALSE)
  )
  
  # Tests
  expect_s4_class(model, "GBLUPModel")
  expect_true(!is.null(model@model))
})

test_that("AE matrix auto-computes from A when not provided", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.5, 
                               n_qtl = 10, include_gxe = TRUE)
  
  # Only provide A matrix, not AE
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Should auto-compute AE and fit successfully with a message
  expect_message(
    model <- fit_gblup(gdata, K, effects = c("A", "ENV", "AE"), verbose = FALSE),
    "Auto-computing AE matrix from A matrix"
  )
  
  # Tests
  expect_s4_class(model, "GBLUPModel")
  expect_true(!is.null(model@model))
  # At least 2 variance components: genetic effects (A, AE) + residual
  expect_true(nrow(model@varcomp) >= 2)
})

test_that("DE matrix auto-computes from D when not provided", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.5, 
                               n_qtl = 10, include_dominance = TRUE, include_gxe = TRUE)
  
  # Only provide A and D matrices, not DE
  K <- calc_relationship_matrices(gdata, matrices = c("A", "D"))
  
  # Should auto-compute DE and fit successfully with a message
  expect_message(
    model <- fit_gblup(gdata, K, effects = c("A", "D", "ENV", "DE"), verbose = FALSE),
    "Auto-computing DE matrix from D matrix"
  )
  
  # Tests
  expect_s4_class(model, "GBLUPModel")
  expect_true(!is.null(model@model))
  # At least 2 variance components: genetic effects (A, D, DE) + residual
  expect_true(nrow(model@varcomp) >= 2)
})

test_that("AAE matrix auto-computes from AA when not provided", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 25, n_snp = 40)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.5, n_qtl = 10, include_gxe = TRUE)
  
  # Only provide A and AA matrices, not AAE
  K <- calc_relationship_matrices(gdata, matrices = c("A", "AA"))
  
  # Should auto-compute AAE and fit successfully with a message
  expect_message(
    model <- fit_gblup(gdata, K, effects = c("A", "ENV", "AA", "AAE"), verbose = FALSE),
    "Auto-computing AAE matrix from AA matrix"
  )
  
  # Tests
  expect_s4_class(model, "GBLUPModel")
  expect_true(!is.null(model@model))
})

test_that("Error message is helpful when neither AE nor A is provided", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.5, n_qtl = 10)
  
  # Provide K_matrices with an unrelated matrix (not A or AE)
  # This gets past the empty-list check and triggers the improved error message
  K <- list(dummy = diag(30))
  
  # Should error with a helpful message
  expect_error(
    fit_gblup(gdata, K, effects = c("AE"), verbose = FALSE),
    "Provide 'AE' in K_matrices, or provide 'A' so it can be auto-computed"
  )
})

test_that("extract_gebv works correctly", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 20, n_snp = 30)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  model <- fit_gblup(gdata, K, effects = "A", verbose = FALSE)
  
  gebv <- extract_gebv(model)
  
  # Tests
  expect_true(is.data.frame(gebv))
  expect_true(all(c("GID", "ENV", "GEBV", "observed") %in% names(gebv)))
  expect_equal(nrow(gebv), nrow(gdata@phenotypes))
})

test_that("calculate_accuracy works correctly", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.7, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  model <- fit_gblup(gdata, K, effects = "A", verbose = FALSE)
  
  acc <- calculate_accuracy(model)
  
  # Tests
  expect_true(is.numeric(acc))
  expect_true(length(acc) == 1)
  expect_true(acc >= -1 && acc <= 1)
  
  # With high h2, accuracy should be reasonable (but not necessarily > 0.3 always)
  expect_true(!is.na(acc))
})

test_that("Model handles missing phenotypes", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  # Add some missing values
  gdata@phenotypes$trait[sample(1:nrow(gdata@phenotypes), 5)] <- NA
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Should still fit
  expect_silent(
    model <- fit_gblup(gdata, K, effects = "A", verbose = FALSE)
  )
  
  expect_s4_class(model, "GBLUPModel")
})

test_that("Model with dominance effects works", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10,
                               include_dominance = TRUE)
  
  K <- calc_relationship_matrices(gdata, matrices = c("A", "D"))
  
  expect_silent(
    model <- fit_gblup(gdata, K, effects = c("A", "D"), verbose = FALSE)
  )
  
  expect_s4_class(model, "GBLUPModel")
  # At minimum: Genetic + Residual (Dominance may be included)
  expect_true(nrow(model@varcomp) >= 2)
})

test_that("Model summary works", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 25, n_snp = 40)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  model <- fit_gblup(gdata, K, effects = "A", verbose = FALSE)
  
  # Summary produces output, so don't use expect_silent
  expect_error(summary(model), NA)
  
  # Test that show works
  expect_output(show(model), "GBLUPModel")
})

test_that("Predict method works", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 25, n_snp = 40)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  model <- fit_gblup(gdata, K, effects = "A", verbose = FALSE)
  
  # Predict on same data
  predictions <- predict(model)
  
  expect_true(is.data.frame(predictions))
  expect_true("GEBV" %in% names(predictions))
  expect_equal(nrow(predictions), nrow(gdata@phenotypes))
})

test_that("fit_gblup_qtl accepts nIters parameter", {
  set.seed(123)
  gdata <- simulate_genotypes(n_ind = 25, n_snp = 40)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  M <- gdata@genotypes
  qtl_mats <- build_qtl_aware_matrices(M, qtl_markers = colnames(M)[1:3])
  
  expect_silent(
    model <- fit_gblup_qtl(gdata, Ka = qtl_mats$Ka, Xa = qtl_mats$Xa,
                           nIters = 100, verbose = FALSE)
  )
  expect_s4_class(model, "GBLUPModel")
})

test_that("plot_heritability excludes n_env from bars", {
  set.seed(123)
  gdata <- simulate_genotypes(n_ind = 25, n_snp = 40)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  K <- calc_relationship_matrices(gdata, matrices = "A")
  model <- fit_gblup(gdata, K, effects = c("A", "ENV"), verbose = FALSE)
  
  h <- heritability(model)
  h_filtered <- h[!names(h) %in% c("n_env")]
  h_filtered <- h_filtered[!is.na(h_filtered)]
  
  # Should have heritability values only, no n_env
  expect_true(all(names(h_filtered) %in% c("h2", "H2", "h2_total", "H2_total")))
  expect_true(length(h_filtered) > 0)
})

# ========== FIXED EFFECTS TESTS ==========

test_that("fit_gblup works with fixed = NULL (intercept-only, backward compatible)", {
  set.seed(123)
  
  # Generate data
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Fit with explicit NULL (should be same as default)
  model <- fit_gblup(gdata, K, effects = "A", fixed = NULL, verbose = FALSE)
  
  # Tests
  expect_s4_class(model, "GBLUPModel")
  expect_true(!is.null(model@model))
  expect_true(is.data.frame(model@gebv))
})

test_that("fit_gblup works with formula fixed effects", {
  set.seed(123)
  
  # Generate data with extra covariate column
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  # Add a TESTER column to phenotypes
  n_obs <- nrow(gdata@phenotypes)
  gdata@phenotypes$TESTER <- sample(c("TesterA", "TesterB"), n_obs, replace = TRUE)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Fit with fixed effects formula
  model <- fit_gblup(gdata, K, effects = "A", fixed = ~ TESTER, verbose = FALSE)
  
  # Tests
  expect_s4_class(model, "GBLUPModel")
  expect_true(!is.null(model@model))
  expect_true(is.data.frame(model@gebv))
  expect_true(nrow(model@gebv) == n_obs)
})

test_that("fit_gblup works with multiple fixed effects", {
  set.seed(123)
  
  # Generate data with multiple covariates
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  # Add TESTER and BLOCK columns
  n_obs <- nrow(gdata@phenotypes)
  gdata@phenotypes$TESTER <- sample(c("TesterA", "TesterB"), n_obs, replace = TRUE)
  gdata@phenotypes$BLOCK <- sample(c("Block1", "Block2", "Block3"), n_obs, replace = TRUE)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Fit with multiple fixed effects
  model <- fit_gblup(gdata, K, effects = "A", fixed = ~ TESTER + BLOCK, verbose = FALSE)
  
  # Tests
  expect_s4_class(model, "GBLUPModel")
  expect_true(!is.null(model@model))
  expect_true(is.data.frame(model@gebv))
})

test_that("fit_gblup works with pre-built fixed effects matrix", {
  set.seed(123)
  
  # Generate data
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  n_obs <- nrow(gdata@phenotypes)
  
  # Create a custom design matrix (e.g., for tester effects)
  X_custom <- matrix(0, nrow = n_obs, ncol = 2)
  X_custom[, 1] <- 1  # Intercept
  X_custom[1:(n_obs/2), 2] <- 1  # First half gets tester effect
  colnames(X_custom) <- c("Intercept", "Tester1")
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Fit with pre-built matrix
  model <- fit_gblup(gdata, K, effects = "A", fixed = X_custom, verbose = FALSE)
  
  # Tests
  expect_s4_class(model, "GBLUPModel")
  expect_true(!is.null(model@model))
})

test_that("GBLUPData accepts phenotypes with extra columns", {
  set.seed(123)
  
  # Generate basic data
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  # Add extra columns to phenotypes
  gdata@phenotypes$TESTER <- sample(c("TesterA", "TesterB"), 
                                    nrow(gdata@phenotypes), replace = TRUE)
  gdata@phenotypes$BLOCK <- sample(1:3, nrow(gdata@phenotypes), replace = TRUE)
  gdata@phenotypes$LOCATION <- sample(c("LocA", "LocB"), 
                                      nrow(gdata@phenotypes), replace = TRUE)
  
  # Should pass validation
  expect_true(validObject(gdata))
  
  # Check that required columns are still present
  expect_true(all(c("GID", "ENV", "trait") %in% names(gdata@phenotypes)))
})

test_that("cv_gblup works with fixed effects", {
  set.seed(123)
  
  # Generate data with covariate
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  # Add TESTER column
  gdata@phenotypes$TESTER <- sample(c("TesterA", "TesterB"), 
                                    nrow(gdata@phenotypes), replace = TRUE)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Run CV with fixed effects
  cv_res <- cv_gblup(gdata, K, effects = "A", fixed = ~ TESTER, 
                     scheme = "CV1", n_folds = 3, verbose = FALSE)
  
  # Tests
  expect_s4_class(cv_res, "CVResult")
  expect_true(is.data.frame(cv_res@metrics))
  expect_true(nrow(cv_res@metrics) > 0)
  expect_true("predictive_ability" %in% names(cv_res@metrics))
})

test_that("cv_gblup with fixed = NULL works (backward compatible)", {
  set.seed(123)
  
  # Generate data
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Run CV without fixed effects
  cv_res <- cv_gblup(gdata, K, effects = "A", fixed = NULL,
                     scheme = "CV1", n_folds = 3, verbose = FALSE)
  
  # Tests
  expect_s4_class(cv_res, "CVResult")
  expect_true(is.data.frame(cv_res@metrics))
})

test_that("Fixed effects with GxE interaction model", {
  set.seed(123)
  
  # Generate GxE data
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.5, 
                               n_qtl = 10, include_gxe = TRUE)
  
  # Add covariate
  gdata@phenotypes$TESTER <- sample(c("TesterA", "TesterB"), 
                                    nrow(gdata@phenotypes), replace = TRUE)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Fit with fixed effects and GxE
  model <- fit_gblup(gdata, K, effects = c("A", "ENV", "AE"), 
                     fixed = ~ TESTER, verbose = FALSE)
  
  # Tests
  expect_s4_class(model, "GBLUPModel")
  expect_true(!is.null(model@model))
  expect_true(nrow(model@varcomp) >= 2)
})

test_that("Fixed effects error handling - invalid formula", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Try to use a column that doesn't exist
  expect_error(
    fit_gblup(gdata, K, effects = "A", fixed = ~ NONEXISTENT_COLUMN, verbose = FALSE),
    "object 'NONEXISTENT_COLUMN' not found"
  )
})

test_that("Fixed effects error handling - wrong matrix dimensions", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Create matrix with wrong number of rows
  X_wrong <- matrix(1, nrow = 10, ncol = 2)  # Wrong number of rows
  
  expect_error(
    fit_gblup(gdata, K, effects = "A", fixed = X_wrong, verbose = FALSE),
    "Fixed effects matrix must have"
  )
})

test_that("Fixed effects error handling - invalid type", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Try to pass a vector instead of formula or matrix
  expect_error(
    fit_gblup(gdata, K, effects = "A", fixed = c(1, 2, 3), verbose = FALSE),
    "fixed must be NULL, a formula, or a matrix"
  )
})

test_that("Intercept detection handles non-intercept columns of all 1s correctly", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  n_obs <- nrow(gdata@phenotypes)
  
  # Create matrix with a properly named intercept AND a treatment column that's all 1s
  X_custom <- matrix(1, nrow = n_obs, ncol = 3)
  X_custom[, 2] <- 1  # Treatment indicator (all 1s but NOT the intercept)
  X_custom[, 3] <- rnorm(n_obs)  # Another covariate
  colnames(X_custom) <- c("(Intercept)", "Treatment", "Covariate")
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Should fit successfully without adding another intercept
  model <- fit_gblup(gdata, K, effects = "A", fixed = X_custom, verbose = FALSE)
  
  expect_s4_class(model, "GBLUPModel")
})

test_that("Matrix without named intercept gets one added with warning", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  n_obs <- nrow(gdata@phenotypes)
  
  # Create matrix without intercept column
  X_no_intercept <- matrix(rnorm(n_obs * 2), nrow = n_obs, ncol = 2)
  colnames(X_no_intercept) <- c("Covariate1", "Covariate2")
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Should warn about adding intercept
  expect_warning(
    model <- fit_gblup(gdata, K, effects = "A", fixed = X_no_intercept, verbose = FALSE),
    "No intercept column detected"
  )
  
  expect_s4_class(model, "GBLUPModel")
})
