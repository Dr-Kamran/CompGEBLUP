# Test GWAS Functions

test_that("GWAS runs correctly", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 50, n_snp = 100)
  gdata <- simulate_phenotypes(gdata, n_env = 1, h2 = 0.6, n_qtl = 5)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  expect_silent(
    gwas_res <- gwas(gdata, K_matrices = K, effects = "A", verbose = FALSE)
  )
  
  # Tests - updated for GWASResult S4 class
  expect_s4_class(gwas_res, "GWASResult")
  expect_true(is.data.frame(gwas_res@results))
  expect_true(all(c("beta", "se", "t_stat", "p_value", "marker") %in% 
                    names(gwas_res@results)))
  
  # Should test markers (or less after MAF filter)
  expect_true(nrow(gwas_res@results) <= 100)
  expect_true(nrow(gwas_res@results) > 0)
  
  # P-values should be between 0 and 1
  expect_true(all(gwas_res@results$p_value >= 0, na.rm = TRUE))
  expect_true(all(gwas_res@results$p_value <= 1, na.rm = TRUE))
})

test_that("GWAS MAF filter works", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 50, n_snp = 100)
  gdata <- simulate_phenotypes(gdata, n_env = 1, h2 = 0.6, n_qtl = 5)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Test with different MAF thresholds
  gwas_res1 <- gwas(gdata, K_matrices = K, min.MAF = 0.01, verbose = FALSE)
  gwas_res2 <- gwas(gdata, K_matrices = K, min.MAF = 0.1, verbose = FALSE)
  
  # Higher MAF threshold should result in fewer SNPs
  expect_true(nrow(gwas_res2@results) <= nrow(gwas_res1@results))
})

test_that("plot_manhattan works", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 40, n_snp = 100)
  gdata <- simulate_phenotypes(gdata, n_env = 1, h2 = 0.6, n_qtl = 5)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  gwas_res <- gwas(gdata, K_matrices = K, verbose = FALSE)
  
  # Should create plot without error
  expect_silent(plot_manhattan(gwas_res))
})

test_that("GWAS without population structure control works", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 40, n_snp = 80)
  gdata <- simulate_phenotypes(gdata, n_env = 1, h2 = 0.6, n_qtl = 5)
  
  # GWAS without K_matrices - will produce message about testing
  gwas_res <- gwas(gdata, K_matrices = NULL, effects = NULL, verbose = FALSE)
  
  expect_s4_class(gwas_res, "GWASResult")
  expect_true(nrow(gwas_res@results) > 0)
})

test_that("GWAS identifies QTL", {
  set.seed(456)
  
  gdata <- simulate_genotypes(n_ind = 60, n_snp = 100)
  # Simulate with known QTL
  gdata <- simulate_phenotypes(gdata, n_env = 1, h2 = 0.8, n_qtl = 3)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  gwas_res <- gwas(gdata, K_matrices = K, verbose = FALSE)
  
  # Should have some significant SNPs
  significant <- sum(gwas_res@results$p_value < 0.01, na.rm = TRUE)
  expect_true(significant > 0)
})

test_that("plot_qq works", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 40, n_snp = 100)
  gdata <- simulate_phenotypes(gdata, n_env = 1, h2 = 0.6, n_qtl = 5)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  gwas_res <- gwas(gdata, K_matrices = K, verbose = FALSE)
  
  # Should create QQ plot without error
  expect_silent(plot_qq(gwas_res))
})

test_that("GWAS handles multi-environment data", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 40, n_snp = 80)
  gdata <- simulate_phenotypes(gdata, n_env = 2, h2 = 0.6, n_qtl = 5)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Should work with multi-env data
  gwas_res <- gwas(gdata, K_matrices = K, effects = c("A", "ENV"), verbose = FALSE)
  
  expect_s4_class(gwas_res, "GWASResult")
})

test_that("get_top_markers works", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 40, n_snp = 100)
  gdata <- simulate_phenotypes(gdata, n_env = 1, h2 = 0.6, n_qtl = 5)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  gwas_res <- gwas(gdata, K_matrices = K, verbose = FALSE)
  
  # Get top 10 markers
  top <- get_top_markers(gwas_res, n = 10)
  
  expect_true(is.data.frame(top))
  expect_true(nrow(top) <= 10)
  expect_true("log10p" %in% names(top))
})

test_that("GWAS show method works", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 40, n_snp = 100)
  gdata <- simulate_phenotypes(gdata, n_env = 1, h2 = 0.6, n_qtl = 5)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  gwas_res <- gwas(gdata, K_matrices = K, verbose = FALSE)
  
  expect_output(show(gwas_res), "GWAS Result")
  expect_output(show(gwas_res), "markers tested")
  expect_output(show(gwas_res), "lambda")
})

test_that("Lambda calculation works", {
  set.seed(123)
  
  gdata <- simulate_genotypes(n_ind = 50, n_snp = 100)
  gdata <- simulate_phenotypes(gdata, n_env = 1, h2 = 0.6, n_qtl = 5)
  
  K <- calc_relationship_matrices(gdata, matrices = "A")
  gwas_res <- gwas(gdata, K_matrices = K, verbose = FALSE)
  
  # Lambda should be numeric and reasonable
  expect_true(is.numeric(gwas_res@lambda))
  expect_true(gwas_res@lambda > 0)
  expect_true(gwas_res@lambda < 5)  # Usually between 0.5 and 2
})