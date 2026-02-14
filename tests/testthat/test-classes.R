test_that("GBLUPData class works correctly", {
  
  M <- matrix(sample(c(-1, 0, 1), 100, replace = TRUE), 10, 10)
  rownames(M) <- paste0("ID", 1:10)
  colnames(M) <- paste0("SNP", 1:10)
  
  pheno <- data.frame(
    GID = rep(paste0("ID", 1:10), 2),
    ENV = rep(c("E1", "E2"), each = 10),
    trait = rnorm(20, 100, 10)
  )
  
  maf <- colMeans(M + 1) / 2
  
  expect_silent(
    gdata <- new("GBLUPData",
                 genotypes = M,
                 phenotypes = pheno,
                 maf = maf)
  )
  
  expect_true(is.matrix(gdata@genotypes))
  expect_true(is.data.frame(gdata@phenotypes))
  expect_equal(nrow(gdata@genotypes), 10)
  expect_output(show(gdata), "GBLUPData")
})

test_that("GBLUPModel class works correctly", {
  M <- matrix(sample(c(-1, 0, 1), 100, replace = TRUE), 10, 10)
  rownames(M) <- paste0("ID", 1:10)
  colnames(M) <- paste0("SNP", 1:10)
  
  pheno <- data.frame(
    GID = paste0("ID", 1:10),
    ENV = "E1",
    trait = rnorm(10, 100, 10)
  )
  
  gdata <- new("GBLUPData",
               genotypes = M,
               phenotypes = pheno,
               maf = colMeans(M + 1) / 2)
  
  gebv <- data.frame(
    GID = paste0("ID", 1:10),
    ENV = "E1",
    GEBV = rnorm(10, 100, 5),
    observed = pheno$trait
  )
  
  varcomp <- data.frame(
    Component = c("A", "Residual"),
    Variance = c(25, 75)
  )
  
  model <- new("GBLUPModel",
               model = list(convergence = TRUE),
               data = gdata,
               gebv = gebv,
               varcomp = varcomp)
  
  expect_s4_class(model, "GBLUPModel")
  expect_output(show(model), "GBLUPModel")
})

test_that("CVResult class works correctly", {
  predictions <- data.frame(
    GID = paste0("ID", 1:50),
    observed = rnorm(50, 100, 10),
    predicted = rnorm(50, 100, 8),
    fold = rep(1:5, each = 10),
    replication = rep(1, 50),
    stringsAsFactors = FALSE
  )
  
  metrics <- data.frame(
    fold = 1:5,
    replication = rep(1, 5),
    predictive_ability = runif(5, 0.3, 0.7),
    mse = runif(5, 50, 150),
    n_test = rep(10, 5),
    stringsAsFactors = FALSE
  )
  
  cv_res <- new("CVResult",
                predictions = predictions,
                metrics = metrics,
                fold_results = list(),
                scheme = "CV1",
                n_folds = 5,
                n_reps = 1,
                effects = "A")
  
  expect_s4_class(cv_res, "CVResult")
  expect_output(show(cv_res), "Cross-Validation")
  expect_output(show(cv_res), "CV1")
})

test_that("GWASResult class works correctly", {
  results <- data.frame(
    marker = paste0("SNP", 1:100),
    beta = rnorm(100, 0, 0.5),
    se = runif(100, 0.1, 0.3),
    t_stat = rnorm(100, 0, 2),
    p_value = runif(100, 0, 1),
    stringsAsFactors = FALSE
  )
  
  gwas_res <- new("GWASResult",
                  results = results,
                  lambda = 1.02,
                  n_markers = 100,
                  min_MAF = 0.05,
                  method = "Wald")
  
  expect_s4_class(gwas_res, "GWASResult")
  expect_output(show(gwas_res), "GWAS Result")
  expect_output(show(gwas_res), "lambda")
})

test_that("Class validation works", {
  # Test GBLUPData validation - missing required columns
  expect_error(
    new("GBLUPData",
        genotypes = matrix(sample(c(-1, 0, 1), 10, replace = TRUE), 5, 2),
        phenotypes = data.frame(x = 1:5),  # Missing GID, ENV, trait
        maf = c(0.3, 0.4)),
    "phenotypes must have columns"
  )
  
  # Test GBLUPData validation - invalid genotype values
  expect_error(
    new("GBLUPData",
        genotypes = matrix(c(0, 1, 2, 3), 2, 2),  # Invalid: should be -1, 0, 1
        phenotypes = data.frame(GID = c("A", "B"), ENV = "E1", trait = c(1, 2)),
        maf = c(0.3, 0.4)),
    "genotypes must be coded as exactly -1, 0, 1"
  )
  
  # Test GBLUPData validation - maf length mismatch
  expect_error(
    new("GBLUPData",
        genotypes = matrix(sample(c(-1, 0, 1), 10, replace = TRUE), 5, 2),
        phenotypes = data.frame(GID = paste0("ID", 1:5), ENV = "E1", trait = rnorm(5)),
        maf = c(0.3, 0.4, 0.5)),  # 3 MAF values but only 2 SNPs
    "maf length must match"
  )
})

test_that("calculate_accuracy generic works", {
  M <- matrix(sample(c(-1, 0, 1), 100, replace = TRUE), 10, 10)
  rownames(M) <- paste0("ID", 1:10)
  colnames(M) <- paste0("SNP", 1:10)
  
  pheno <- data.frame(
    GID = paste0("ID", 1:10),
    ENV = "E1",
    trait = rnorm(10, 100, 10)
  )
  
  gdata <- new("GBLUPData",
               genotypes = M,
               phenotypes = pheno,
               maf = colMeans(M + 1) / 2)
  
  # Create correlated GEBV and observed values
  observed <- rnorm(10, 100, 10)
  gebv_vals <- observed * 0.8 + rnorm(10, 0, 3)  # Correlated
  
  gebv <- data.frame(
    GID = paste0("ID", 1:10),
    ENV = "E1",
    GEBV = gebv_vals,
    observed = observed
  )
  
  varcomp <- data.frame(
    Component = c("A", "Residual"),
    Variance = c(25, 75)
  )
  
  model <- new("GBLUPModel",
               model = list(convergence = TRUE),
               data = gdata,
               gebv = gebv,
               varcomp = varcomp)
  
  acc <- calculate_accuracy(model)
  
  expect_true(is.numeric(acc))
  expect_true(length(acc) == 1)
  expect_true(acc >= -1 && acc <= 1)
})

test_that("extract_gebv and extract_varcomp work", {
  M <- matrix(sample(c(-1, 0, 1), 100, replace = TRUE), 10, 10)
  rownames(M) <- paste0("ID", 1:10)
  colnames(M) <- paste0("SNP", 1:10)
  
  pheno <- data.frame(
    GID = paste0("ID", 1:10),
    ENV = "E1",
    trait = rnorm(10, 100, 10)
  )
  
  gdata <- new("GBLUPData",
               genotypes = M,
               phenotypes = pheno,
               maf = colMeans(M + 1) / 2)
  
  gebv <- data.frame(
    GID = paste0("ID", 1:10),
    ENV = "E1",
    GEBV = rnorm(10, 100, 5),
    observed = pheno$trait
  )
  
  varcomp <- data.frame(
    Component = c("A", "Residual"),
    Variance = c(25, 75)
  )
  
  model <- new("GBLUPModel",
               model = list(convergence = TRUE),
               data = gdata,
               gebv = gebv,
               varcomp = varcomp)
  
  # Test extract_gebv
  gebv_extracted <- extract_gebv(model)
  expect_true(is.data.frame(gebv_extracted))
  expect_true(all(c("GID", "ENV", "GEBV", "observed") %in% names(gebv_extracted)))
  
  # Test extract_varcomp
  vc_extracted <- extract_varcomp(model)
  expect_true(is.data.frame(vc_extracted))
  expect_true(all(c("Component", "Variance") %in% names(vc_extracted)))
})