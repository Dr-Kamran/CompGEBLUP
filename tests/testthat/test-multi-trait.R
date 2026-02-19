test_that("Multi-trait model fits successfully", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 50, n_snp = 100)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  
  # Add second trait
  gdata@phenotypes$trait2 <- gdata@phenotypes$trait * 0.8 + rnorm(nrow(gdata@phenotypes), 0, 5)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  mt_model <- fit_multi_trait(
    gdata,
    traits = c("trait", "trait2"),
    K_matrices = K,
    effects = "A",
    genetic_model = "diagonal",
    verbose = FALSE
  )
  
  # Tests
  expect_true(!is.null(mt_model$model))
  expect_true(is.matrix(mt_model$genetic_correlation))
  expect_equal(nrow(mt_model$genetic_correlation), 2)
  expect_equal(length(mt_model$heritability), 2)
})

test_that("Multi-trait genetic correlations are bounded", {
  set.seed(456)
  
  gdata <- simulate_genotypes(n_ind = 50, n_snp = 100)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  gdata@phenotypes$trait2 <- gdata@phenotypes$trait + rnorm(nrow(gdata@phenotypes), 0, 10)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  mt_model <- fit_multi_trait(
    gdata,
    traits = c("trait", "trait2"),
    K_matrices = K,
    effects = "A",
    genetic_model = "unstructured",  # Full correlation matrix
    verbose = FALSE
  )
  
  # Genetic correlations should be between -1 and 1
  gen_cor <- mt_model$genetic_correlation
  off_diag <- gen_cor[lower.tri(gen_cor) | upper.tri(gen_cor)]
  
  expect_true(all(abs(off_diag) <= 1, na.rm = TRUE))
})

test_that("Selection index calculation works", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 40, n_snp = 80)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 10)
  gdata@phenotypes$trait2 <- gdata@phenotypes$trait * 0.7 + rnorm(nrow(gdata@phenotypes), 0, 5)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  mt_model <- fit_multi_trait(gdata, c("trait", "trait2"), K, "A", verbose = FALSE)
  
  weights <- c(trait = 0.6, trait2 = 0.4)
  index <- selection_index(mt_model, weights)
  
  expect_true(is.data.frame(index))
  expect_true("Index" %in% names(index))
  expect_equal(nrow(index), length(unique(gdata@phenotypes$GID)))
})