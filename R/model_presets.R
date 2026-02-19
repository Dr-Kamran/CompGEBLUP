#' Fit Pre-defined GBLUP Models
#'
#' Convenience functions to fit the 5 core GBLUP models.
#' Updated in v2.3.1 to reflect empirically validated model structures
#' from G2F 2021 analysis:
#' \itemize{
#'   \item Model 4 drops DE by default (DE causes severe overfitting in hybrid data)
#'   \item Model 5 drops D by default (σ²_D often hits boundary in multi-effect models)
#'   \item All presets use \code{ridge_per_matrix} for differential regularisation
#'   \item Default \code{nIters = 500} (50 is too few for multi-effect convergence)
#' }
#'
#' @param gdata GBLUPData object
#' @param ridge Base ridge parameter (default: 0.001). Used for well-conditioned matrices.
#' @param ridge_D Ridge for dominance matrix (default: 0.05). D has many near-zero eigenvalues.
#' @param ridge_AA Ridge for epistasis matrix (default: 0.05). AA Hadamard is ill-conditioned.
#' @param ridge_AE Ridge for GxE Kronecker matrix (default: 0.05). AE is high-dimensional.
#' @param nIters Maximum REML iterations (default: 500). Increase for complex models.
#' @param verbose Logical, print progress messages
#' @param fixed Optional fixed effects specification (e.g., ~ TESTER). Default is NULL.
#' @param use.emma Solver selection (default: FALSE for multi-effect, TRUE for M1).
#' @param include_DE Logical, include dominance-by-environment interaction in Model 4
#'   (default: FALSE). DE typically causes overfitting in hybrid breeding data.
#' @param include_D Logical, include dominance in Model 5 (default: FALSE).
#'   σ²_D is often a boundary estimate when epistasis is already in the model.
#' @param include_AAE Logical, include epistasis-by-environment in Model 5 (default: FALSE).
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
#' m4 <- fit_model4_ADE(gdata)               # No DE by default
#' m4_de <- fit_model4_ADE(gdata, include_DE = TRUE)  # With DE
#' m5 <- fit_model5_ADAA(gdata)              # No D by default
#'
#' # Compare models
#' compare_models(list(M1 = m1, M2 = m2, M3 = m3, M4 = m4, M5 = m5))
#' }
NULL

#' @rdname model_presets
#' @export
fit_model1_A <- function(gdata, ridge = 0.001, nIters = 500, verbose = TRUE,
                         fixed = NULL) {
  if (verbose) cat("=== Fitting Model 1: A ===\n")
  K <- calc_relationship_matrices(gdata, matrices = "A", ridge = ridge)
  model <- fit_gblup(gdata, K, effects = "A", ridge = ridge,
                     nIters = nIters, verbose = verbose, fixed = fixed,
                     use.emma = TRUE)  # Single effect: EMMA is optimal
  model@description <- "Model 1: Additive (A)"
  return(model)
}

#' @rdname model_presets
#' @export
fit_model2_AD <- function(gdata, ridge = 0.001, ridge_D = 0.05,
                          nIters = 500, verbose = TRUE, fixed = NULL) {
  if (verbose) cat("=== Fitting Model 2: A + D ===\n")
  K <- calc_relationship_matrices(gdata, matrices = c("A", "D"), ridge = ridge)
  model <- fit_gblup(gdata, K, effects = c("A", "D"),
                     ridge = ridge,
                     ridge_per_matrix = c(A = ridge, D = ridge_D),
                     nIters = nIters, verbose = verbose, fixed = fixed,
                     use.emma = FALSE)  # Multi-effect: Henderson EM-REML
  model@description <- "Model 2: Additive + Dominance (A + D)"
  return(model)
}

#' @rdname model_presets
#' @export
fit_model3_AE <- function(gdata, ridge = 0.001, ridge_AE = 0.05,
                          nIters = 500, verbose = TRUE, fixed = NULL) {
  if (verbose) cat("=== Fitting Model 3: A + ENV + AE ===\n")
  K <- calc_relationship_matrices(gdata, matrices = c("A", "AE"), ridge = ridge)
  model <- fit_gblup(gdata, K, effects = c("A", "ENV", "AE"),
                     ridge = ridge,
                     ridge_per_matrix = c(A = ridge, AE = ridge_AE),
                     nIters = nIters, verbose = verbose, fixed = fixed,
                     use.emma = FALSE)
  model@description <- "Model 3: Additive GxE (A + ENV + AE)"
  return(model)
}

#' @rdname model_presets
#' @export
fit_model4_ADE <- function(gdata, include_DE = FALSE,
                           ridge = 0.001, ridge_D = 0.05, ridge_AE = 0.05,
                           nIters = 500, verbose = TRUE, fixed = NULL) {
  if (include_DE) {
    if (verbose) cat("=== Fitting Model 4: A + D + ENV + AE + DE ===\n")
    K <- calc_relationship_matrices(gdata, matrices = c("A", "D", "AE", "DE"),
                                    ridge = ridge)
    eff <- c("A", "D", "ENV", "AE", "DE")
    rpm <- c(A = ridge, D = ridge_D, AE = ridge_AE, DE = ridge_AE)
    desc <- "Model 4: Full GxE with Dominance (A + D + ENV + AE + DE)"
  } else {
    if (verbose) cat("=== Fitting Model 4: A + D + ENV + AE (DE dropped) ===\n")
    K <- calc_relationship_matrices(gdata, matrices = c("A", "D", "AE"),
                                    ridge = ridge)
    eff <- c("A", "D", "ENV", "AE")
    rpm <- c(A = ridge, D = ridge_D, AE = ridge_AE)
    desc <- "Model 4: GxE with Dominance, no DE (A + D + ENV + AE)"
  }
  
  model <- fit_gblup(gdata, K, effects = eff,
                     ridge = ridge, ridge_per_matrix = rpm,
                     nIters = nIters, verbose = verbose, fixed = fixed,
                     use.emma = FALSE)
  model@description <- desc
  return(model)
}

#' @rdname model_presets
#' @export
fit_model5_ADAA <- function(gdata, include_D = FALSE, include_AAE = FALSE,
                            ridge = 0.001, ridge_D = 0.05, ridge_AA = 0.05,
                            ridge_AE = 0.05, nIters = 500, verbose = TRUE,
                            fixed = NULL) {
  # Build effect list and matrices dynamically
  mat_names <- "A"
  eff <- "A"
  rpm <- c(A = ridge)
  
  if (include_D) {
    mat_names <- c(mat_names, "D")
    eff <- c(eff, "D")
    rpm <- c(rpm, D = ridge_D)
  }
  
  # Always include ENV + AE for GxE
  mat_names <- c(mat_names, "AE")
  eff <- c(eff, "ENV", "AE")
  rpm <- c(rpm, AE = ridge_AE)
  
  # Always include AA (epistasis)
  mat_names <- c(mat_names, "AA")
  eff <- c(eff, "AA")
  rpm <- c(rpm, AA = ridge_AA)
  
  if (include_AAE) {
    mat_names <- c(mat_names, "AAE")
    eff <- c(eff, "AAE")
    rpm <- c(rpm, AAE = ridge_AE)
  }
  
  if (verbose) cat("=== Fitting Model 5:", paste(eff, collapse = " + "), "===\n")
  
  K <- calc_relationship_matrices(gdata, matrices = mat_names, ridge = ridge)
  model <- fit_gblup(gdata, K, effects = eff,
                     ridge = ridge, ridge_per_matrix = rpm,
                     nIters = nIters, verbose = verbose, fixed = fixed,
                     use.emma = FALSE)
  model@description <- paste("Model 5:", paste(eff, collapse = " + "))
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