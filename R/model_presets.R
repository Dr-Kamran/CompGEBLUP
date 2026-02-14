#' Fit Pre-defined GBLUP Models
#'
#' Convenience functions to fit the 5 core GBLUP models
#'
#' @param gdata GBLUPData object
#' @param include_AAE Logical, include epistasis-by-environment interactions (Model 5 only, default: FALSE)
#' @param ridge Ridge parameter for numerical stability (default: 0.001)
#' @param nIters Maximum REML iterations (default: 50)
#' @param verbose Logical, print progress messages
#' @param fixed Optional fixed effects specification (e.g., ~ TESTER). Default is NULL (intercept only).
#'
#' @return GBLUPModel object
#' @name model_presets
#' @examples
#' \dontrun{
#' gdata <- simulate_genotypes(n_ind = 100, n_snp = 500)
#' gdata <- simulate_phenotypes(gdata, n_env = 3, h2 = 0.6, n_qtl = 20)
#'
#' # Fit each of the 5 core models
#' m1 <- fit_model1_A(gdata)
#' m2 <- fit_model2_AD(gdata)
#' m3 <- fit_model3_AE(gdata)
#' m4 <- fit_model4_ADE(gdata)
#' m5 <- fit_model5_ADAA(gdata)
#'
#' # Compare models
#' compare_models(list(M1 = m1, M2 = m2, M3 = m3, M4 = m4, M5 = m5))
#' }
NULL

#' @rdname model_presets
#' @export
fit_model1_A <- function(gdata, ridge = 0.001, nIters = 50, verbose = TRUE, fixed = NULL) {
  if (verbose) cat("=== Fitting Model 1: A ===\n")
  K <- calc_relationship_matrices(gdata, matrices = "A", ridge = ridge)
  model <- fit_gblup(gdata, K, effects = "A", ridge = ridge, nIters = nIters, verbose = verbose, fixed = fixed)
  model@description <- "Model 1: Additive (A)"
  return(model)
}

#' @rdname model_presets
#' @export
fit_model2_AD <- function(gdata, ridge = 0.001, nIters = 50, verbose = TRUE, fixed = NULL) {
  if (verbose) cat("=== Fitting Model 2: A + D ===\n")
  K <- calc_relationship_matrices(gdata, matrices = c("A", "D"), ridge = ridge)
  model <- fit_gblup(gdata, K, effects = c("A", "D"), ridge = ridge, nIters = nIters, verbose = verbose, fixed = fixed)
  model@description <- "Model 2: Additive + Dominance (A + D)"
  return(model)
}

#' @rdname model_presets
#' @export
fit_model3_AE <- function(gdata, ridge = 0.001, nIters = 50, verbose = TRUE, fixed = NULL) {
  if (verbose) cat("=== Fitting Model 3: A + ENV + AE ===\n")
  K <- calc_relationship_matrices(gdata, matrices = c("A", "AE"), ridge = ridge)
  model <- fit_gblup(gdata, K, effects = c("A", "ENV", "AE"), ridge = ridge, nIters = nIters, verbose = verbose, fixed = fixed)
  model@description <- "Model 3: Additive GxE (A + ENV + AE)"
  return(model)
}

#' @rdname model_presets
#' @export
fit_model4_ADE <- function(gdata, ridge = 0.001, nIters = 50, verbose = TRUE, fixed = NULL) {
  if (verbose) cat("=== Fitting Model 4: A + D + ENV + AE + DE ===\n")
  K <- calc_relationship_matrices(gdata, matrices = c("A", "D", "AE", "DE"), ridge = ridge)
  model <- fit_gblup(gdata, K, effects = c("A", "D", "ENV", "AE", "DE"), ridge = ridge, nIters = nIters, verbose = verbose, fixed = fixed)
  model@description <- "Model 4: Full GxE with Dominance (A + D + ENV + AE + DE)"
  return(model)
}

#' @rdname model_presets
#' @export
fit_model5_ADAA <- function(gdata, include_AAE = FALSE, ridge = 0.001, nIters = 50, verbose = TRUE, fixed = NULL) {
  if (verbose) {
    if (include_AAE) {
      cat("=== Fitting Model 5: A + D + ENV + AE + AA + AAE ===\n")
    } else {
      cat("=== Fitting Model 5: A + D + ENV + AE + AA ===\n")
    }
  }
  
  if (include_AAE) {
    K <- calc_relationship_matrices(gdata, matrices = c("A", "D", "AE", "AA", "AAE"), ridge = ridge)
    model <- fit_gblup(gdata, K, effects = c("A", "D", "ENV", "AE", "AA", "AAE"), 
                       ridge = ridge, nIters = nIters, verbose = verbose, fixed = fixed)
    model@description <- "Model 5: Full GxE with Epistasis (A + D + ENV + AE + AA + AAE)"
  } else {
    K <- calc_relationship_matrices(gdata, matrices = c("A", "D", "AE", "AA"), ridge = ridge)
    model <- fit_gblup(gdata, K, effects = c("A", "D", "ENV", "AE", "AA"), 
                       ridge = ridge, nIters = nIters, verbose = verbose, fixed = fixed)
    model@description <- "Model 5: Full GxE with Epistasis (A + D + ENV + AE + AA)"
  }
  
  return(model)
}

#' Compare Multiple GBLUP Models
#'
#' Compare prediction accuracy and variance components across models
#'
#' @param models Named list of GBLUPModel objects
#' @return Data frame with model comparison statistics
#' @export
#' @examples
#' \dontrun{
#' m1 <- fit_model1_A(gdata, verbose = FALSE)
#' m2 <- fit_model2_AD(gdata, verbose = FALSE)
#' m3 <- fit_model3_AE(gdata, verbose = FALSE)
#'
#' comparison <- compare_models(list(
#'   "Model 1 (A)" = m1,
#'   "Model 2 (A+D)" = m2,
#'   "Model 3 (A+ENV+AE)" = m3
#' ))
#' print(comparison)
#' }
compare_models <- function(models) {
  
  if (!is.list(models) || is.null(names(models))) {
    stop("models must be a named list of GBLUPModel objects")
  }
  
  # Validate all are GBLUPModel objects
  for (nm in names(models)) {
    if (!inherits(models[[nm]], "GBLUPModel")) {
      stop("'", nm, "' is not a GBLUPModel object")
    }
  }
  
  results <- data.frame(
    Model = names(models),
    Predictive_Ability = sapply(models, function(m) {
      tryCatch(calculate_accuracy(m), error = function(e) NA)
    }),
    n_VarianceComponents = sapply(models, function(m) nrow(m@varcomp)),
    Total_Variance = sapply(models, function(m) sum(m@varcomp$Variance)),
    Residual_Variance = sapply(models, function(m) {
      res_idx <- grep("Residual", m@varcomp$Component, ignore.case = TRUE)
      if (length(res_idx) > 0) m@varcomp$Variance[res_idx] else NA
    }),
    LogLik = sapply(models, function(m) {
      if (!is.null(m@model$llik)) m@model$llik else NA
    }),
    AIC = sapply(models, function(m) {
      if (!is.null(m@model$AIC)) m@model$AIC else NA
    }),
    BIC = sapply(models, function(m) {
      if (!is.null(m@model$BIC)) m@model$BIC else NA
    }),
    Converged = sapply(models, function(m) {
      if (!is.null(m@model$convergence)) m@model$convergence else NA
    }),
    stringsAsFactors = FALSE
  )
  
  # Calculate proportion of variance explained
  results$PropVarianceExplained <- 1 - (results$Residual_Variance / results$Total_Variance)
  
  # Rank by predictive ability
  results <- results[order(-results$Predictive_Ability), ]
  rownames(results) <- NULL
  
  return(results)
}