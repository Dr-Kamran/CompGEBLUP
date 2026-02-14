# Test Heritability Function

test_that("ENV variance is not classified as genetic variance", {
  set.seed(123)
  
  # Generate data with multiple environments
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 3, h2 = 0.6, n_qtl = 10)
  
  # Calculate matrices
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Fit model without ENV effect
  model_no_env <- fit_gblup(gdata, K, effects = "A", verbose = FALSE)
  h2_no_env <- heritability(model_no_env)
  
  # Fit model with ENV effect
  model_with_env <- fit_gblup(gdata, K, effects = c("A", "ENV"), verbose = FALSE)
  h2_with_env <- heritability(model_with_env)
  
  # Test 1: Both heritability values should be reasonable (between 0 and 1)
  expect_true(h2_no_env["h2"] > 0 && h2_no_env["h2"] <= 1)
  expect_true(h2_with_env["h2"] > 0 && h2_with_env["h2"] <= 1)
  
  # Test 2: Heritability with ENV should NOT be dramatically higher than without ENV
  # (it should be lower or similar, not inflated)
  # Allow some variation due to estimation, but not more than 20% increase
  expect_true(h2_with_env["h2"] <= h2_no_env["h2"] * 1.2)
  
  # Test 3: ENV variance should be in variance components when included
  expect_true("ENV" %in% model_with_env@varcomp$Component)
  
  # Test 4: ENV variance should be positive
  env_var <- model_with_env@varcomp$Variance[model_with_env@varcomp$Component == "ENV"]
  expect_true(env_var > 0)
})

test_that("n_env is returned as attribute, not in vector", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 3, h2 = 0.6, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  model <- fit_gblup(gdata, K, effects = c("A", "ENV"), verbose = FALSE)
  
  h2_result <- heritability(model)
  
  # Test 1: Result should have exactly 4 elements
  expect_equal(length(h2_result), 4)
  
  # Test 2: Names should be h2, H2, h2_total, H2_total only
  expect_equal(names(h2_result), c("h2", "H2", "h2_total", "H2_total"))
  
  # Test 3: n_env should NOT be in the names
  expect_false("n_env" %in% names(h2_result))
  
  # Test 4: n_env should be accessible via attr()
  n_env_attr <- attr(h2_result, "n_env")
  expect_false(is.null(n_env_attr))
  expect_equal(n_env_attr, 3)
  
  # Test 5: All heritability values should be numeric (not n_env integer)
  expect_true(all(sapply(h2_result, is.numeric)))
  expect_true(all(h2_result >= 0 & h2_result <= 1, na.rm = TRUE))
})

test_that("Heritability calculation handles edge cases correctly", {
  set.seed(456)
  
  # Test with single environment
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 1, h2 = 0.5, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  model <- fit_gblup(gdata, K, effects = "A", verbose = FALSE)
  
  h2_result <- heritability(model)
  
  # Test 1: Should work with single environment
  expect_false(is.na(h2_result["h2"]))
  expect_false(is.na(h2_result["H2"]))
  
  # Test 2: n_env should be 1
  expect_equal(attr(h2_result, "n_env"), 1)
  
  # Test 3: Heritability should be reasonable
  expect_true(h2_result["h2"] > 0 && h2_result["h2"] <= 1)
})

test_that("Kronecker consistency: auto-computed AE matches pre-computed", {
  set.seed(789)
  
  # Generate data
  gdata <- simulate_genotypes(n_ind = 20, n_snp = 40)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.5, 
                               n_qtl = 10, include_gxe = TRUE)
  
  # Fit model with pre-computed AE matrix
  K_full <- calc_relationship_matrices(gdata, matrices = c("A", "AE"))
  model_precomputed <- fit_gblup(gdata, K_full, effects = c("A", "ENV", "AE"), 
                                 verbose = FALSE)
  
  # Fit model with auto-computed AE matrix (only provide A)
  K_partial <- calc_relationship_matrices(gdata, matrices = "A")
  suppressMessages(
    model_autocomputed <- fit_gblup(gdata, K_partial, effects = c("A", "ENV", "AE"), 
                                    verbose = FALSE)
  )
  
  # Extract heritability estimates
  h2_pre <- heritability(model_precomputed)
  h2_auto <- heritability(model_autocomputed)
  
  # Test 1: Both models should fit successfully
  expect_s4_class(model_precomputed, "GBLUPModel")
  expect_s4_class(model_autocomputed, "GBLUPModel")
  
  # Test 2: Heritability estimates should be similar (within 10%)
  # This tests that the Kronecker product order is consistent
  expect_true(abs(h2_pre["h2"] - h2_auto["h2"]) / h2_pre["h2"] < 0.1)
  expect_true(abs(h2_pre["H2"] - h2_auto["H2"]) / h2_pre["H2"] < 0.1)
  
  # Test 3: GEBVs should be correlated
  gebv_pre <- model_precomputed@gebv$GEBV
  gebv_auto <- model_autocomputed@gebv$GEBV
  correlation <- cor(gebv_pre, gebv_auto, use = "complete.obs")
  expect_true(correlation > 0.95)  # Should be very similar
})
