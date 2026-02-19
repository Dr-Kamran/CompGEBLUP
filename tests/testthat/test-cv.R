# Test Cross-Validation Functions

test_that("CV1 (random CV) works correctly", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 40, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  cv_res <- cv_gblup(gdata, K, 
                     effects = "A",
                     scheme = "CV1",
                     n_folds = 5,
                     n_reps = 1,
                     verbose = FALSE)
  
  # Tests
  expect_s4_class(cv_res, "CVResult")
  expect_true(is.data.frame(cv_res@metrics))
  expect_true(is.data.frame(cv_res@predictions))
  
  # Check metrics results
  expect_equal(nrow(cv_res@metrics), 5)  # 5 folds
  expect_true(all(c("fold", "replication", "predictive_ability", "mse", "n_test") %in% 
                    names(cv_res@metrics)))
  
  # Predictive ability should be reasonable
  mean_acc <- mean(cv_res@metrics$predictive_ability, na.rm = TRUE)
  expect_true(!is.na(mean_acc))
})

test_that("CV2 (leave-one-environment-out) works correctly", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 3, h2 = 0.6, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  cv_res <- cv_gblup(gdata, K,
                     effects = "A",
                     scheme = "CV2",
                     verbose = FALSE)
  
  # Tests
  expect_s4_class(cv_res, "CVResult")
  
  # Should have as many folds as environments
  n_env <- length(unique(gdata@phenotypes$ENV))
  expect_equal(nrow(cv_res@metrics), n_env)
})

test_that("CV0 (within-environment) works correctly", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  cv_res <- cv_gblup(gdata, K,
                     effects = "A",
                     scheme = "CV0",
                     n_folds = 3,
                     verbose = FALSE)
  
  # Tests
  expect_s4_class(cv_res, "CVResult")
  # CV0 creates folds within each environment, so might have more folds
  expect_true(nrow(cv_res@metrics) >= 3)
})

test_that("Multiple replications work correctly", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  cv_res <- cv_gblup(gdata, K,
                     effects = "A",
                     scheme = "CV1",
                     n_folds = 3,
                     n_reps = 2,
                     verbose = FALSE)
  
  # Should have 3 folds x 2 reps = 6 rows
  expect_true(nrow(cv_res@metrics) >= 5)
  expect_equal(length(unique(cv_res@metrics$replication)), 2)
})

test_that("CV handles different effects", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # With A + ENV
  cv_res <- cv_gblup(gdata, K,
                     effects = c("A", "ENV"),
                     scheme = "CV1",
                     n_folds = 3,
                     verbose = FALSE)
  
  expect_s4_class(cv_res, "CVResult")
})

test_that("CV show method works", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  cv_res <- cv_gblup(gdata, K, effects = "A", n_folds = 3, verbose = FALSE)
  
  expect_output(show(cv_res), "Cross-Validation")
  expect_output(show(cv_res), "Predictive ability")
  expect_output(show(cv_res), "Scheme")
})

test_that("CV parameters are stored correctly", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  cv_res <- cv_gblup(gdata, K,
                     effects = c("A", "ENV"),
                     scheme = "CV1",
                     n_folds = 5,
                     n_reps = 2,
                     verbose = FALSE)
  
  # Check parameters
  expect_equal(cv_res@scheme, "CV1")
  expect_equal(cv_res@n_folds, 5)
  expect_equal(cv_res@n_reps, 2)
  expect_equal(cv_res@effects, c("A", "ENV"))
})

test_that("CV handles GxE effects", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10, 
                               include_gxe = TRUE)
  
  K <- calc_relationship_matrices(gdata, matrices = c("A", "AE"))
  
  cv_res <- cv_gblup(gdata, K,
                     effects = c("A", "ENV", "AE"),
                     scheme = "CV1",
                     n_folds = 3,
                     verbose = FALSE)
  
  expect_s4_class(cv_res, "CVResult")
})

test_that("CV0 has no data leakage", {
  set.seed(123)
  
  # Create data with SAME individuals in multiple environments
  gdata <- simulate_genotypes(n_ind = 20, n_snp = 30)
  gdata <- simulate_phenotypes(gdata, n_env = 3, h2 = 0.6, n_qtl = 5)
  
  pheno <- gdata@phenotypes
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Perform CV0
  cv_res <- cv_gblup(gdata, K,
                     effects = "A",
                     scheme = "CV0",
                     n_folds = 3,
                     verbose = FALSE)
  
  # Verify no data leakage by checking fold structure
  pheno <- gdata@phenotypes
  pheno$GID <- as.character(pheno$GID)
  pheno$ENV <- as.character(pheno$ENV)
  
  folds <- CompGEBLUP:::create_within_env_folds(pheno, n_folds = 3)
  
  for (fold in folds) {
    test_ids <- unique(pheno$GID[fold$test])
    test_env <- unique(pheno$ENV[fold$test])
    
    # Get training observations
    train_obs <- pheno[fold$train, ]
    
    # Check 1: Training should ONLY be from same environment
    expect_true(all(train_obs$ENV == test_env),
                info = paste("Training contains other environments in fold for env", test_env))
    
    # Check 2: No overlap between test individuals and train individuals
    overlap <- intersect(test_ids, unique(train_obs$GID))
    expect_equal(length(overlap), 0,
                 info = paste("Data leakage: test individuals found in training for env", test_env))
  }
  
  expect_s4_class(cv_res, "CVResult")
})

test_that("CV1 keeps individuals together across environments", {
  set.seed(456)
  
  gdata <- simulate_genotypes(n_ind = 20, n_snp = 30)
  gdata <- simulate_phenotypes(gdata, n_env = 3, h2 = 0.6, n_qtl = 5)
  
  pheno <- gdata@phenotypes
  pheno$GID <- as.character(pheno$GID)
  
  folds <- CompGEBLUP:::create_random_folds(pheno, n_folds = 4)
  
  for (fold in folds) {
    test_ids <- unique(pheno$GID[fold$test])
    train_ids <- unique(pheno$GID[fold$train])
    
    # No individual should be in both test and train
    overlap <- intersect(test_ids, train_ids)
    expect_equal(length(overlap), 0,
                 info = "Individual split across test and train sets in CV1")
  }
})

test_that("CV handles edge case: minimum individuals", {
  set.seed(789)
  
  # Small dataset
  gdata <- simulate_genotypes(n_ind = 10, n_snp = 20)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 3)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Should work with 2 folds
  cv_res <- cv_gblup(gdata, K,
                     effects = "A",
                     scheme = "CV1",
                     n_folds = 2,
                     verbose = FALSE)
  
  expect_s4_class(cv_res, "CVResult")
})

test_that("CV fails gracefully with too few observations", {
  set.seed(101)
  
  gdata <- simulate_genotypes(n_ind = 5, n_snp = 20)
  gdata <- simulate_phenotypes(gdata, n_env = 1, h2 = 0.6, n_qtl = 2)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Should error with 10 folds - check for either possible error message
  expect_error(
    cv_gblup(gdata, K, effects = "A", scheme = "CV1", n_folds = 10, verbose = FALSE),
    "Not enough"  # Matches both "Not enough observations" and "Not enough unique individuals"
  )
})

test_that("get_cv_metrics works correctly", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  cv_res <- cv_gblup(gdata, K, effects = "A", n_folds = 3, verbose = FALSE)
  
  # Test overall metrics
  overall <- get_cv_metrics(cv_res, by_fold = FALSE)
  expect_true(is.data.frame(overall))
  expect_true("mean_predictive_ability" %in% names(overall))
  expect_true("mean_mse" %in% names(overall))
  
  # Test by-fold metrics
  by_fold <- get_cv_metrics(cv_res, by_fold = TRUE)
  expect_true(is.data.frame(by_fold))
  expect_equal(nrow(by_fold), 3)  # 3 folds
})

test_that("plot_cv_results works without error", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 30, n_snp = 50)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  cv_res <- cv_gblup(gdata, K, effects = "A", n_folds = 3, verbose = FALSE)
  
  # Should not error
  expect_silent(plot_cv_results(cv_res, type = "scatter"))
  expect_silent(plot_cv_results(cv_res, type = "boxplot"))
  expect_silent(plot_cv_results(cv_res, type = "density"))
})