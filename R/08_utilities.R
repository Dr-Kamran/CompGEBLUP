# =============================================================================
# Error Handling Patterns for CompGEBLUP
# =============================================================================
# 
# Consistent error handling conventions:
# 1. Input validation: stop() with informative message
# 2. Numerical issues: warning() + return degraded result (with attr explaining why)
# 3. Missing data: Return NA with attr("reason") explaining why
# 4. Recoverable issues: message() for information only
#
# Helper functions below implement these patterns.
# =============================================================================

#' Create Result with Warning Attribute
#' 
#' Returns a value with an attached warning attribute for tracking issues.
#' @param value The result value
#' @param reason Character string explaining the issue
#' @return value with attr("warning_reason")
#' @keywords internal
result_with_warning <- function(value, reason) {
  attr(value, "warning_reason") <- reason
  value
}

#' Create NA Result with Reason
#'
#' Returns NA with an attribute explaining why.
#' @param reason Character string explaining why NA was returned
#' @return NA with attr("na_reason")
#' @keywords internal
na_with_reason <- function(reason) {
  result <- NA
  attr(result, "na_reason") <- reason
  result
}

#' Add Ridge to Matrix Diagonal
#'
#' Adds a small value to the diagonal of a matrix for numerical stability
#'
#' @param K Square matrix
#' @param ridge Ridge value to add to diagonal
#'
#' @return Matrix with ridge added to diagonal
#' @export
#' @examples
#' \dontrun{
#' K <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' K_ridge <- add_ridge(K, ridge = 0.01)
#' }
add_ridge <- function(K, ridge = 0.001) {
  if (!is.matrix(K)) {
    stop("K must be a matrix")
  }
  if (nrow(K) != ncol(K)) {
    stop("K must be a square matrix")
  }
  stopifnot("ridge must be non-negative" = ridge >= 0)
  
  K + diag(ridge, nrow(K))
}

#' Summary Statistics for GBLUPData
#'
#' Prints comprehensive summary of a GBLUPData object.
#'
#' @param gblup_data GBLUPData object
#' @return Invisible NULL
#' @export
#' @examples
#' \dontrun{
#' gdata <- simulate_genotypes(n_ind = 100, n_snp = 500)
#' gdata <- simulate_phenotypes(gdata, n_env = 3, h2 = 0.6)
#' summary_gblup_data(gdata)
#' }
summary_gblup_data <- function(gblup_data) {
  
  if (!inherits(gblup_data, "GBLUPData")) {
    stop("gblup_data must be a GBLUPData object")
  }
  
  cat("GBLUPData Summary\n")
  cat("=================\n\n")
  
  # Genotype summary
  cat("Genotypes:\n")
  cat("  Individuals:", nrow(gblup_data@genotypes), "\n")
  cat("  SNPs:", ncol(gblup_data@genotypes), "\n")
  cat("  MAF range:", paste(round(range(gblup_data@maf), 3), collapse = " - "), "\n")
  cat("  Missing rate:", 
      round(mean(is.na(gblup_data@genotypes)) * 100, 2), "%\n\n")
  
  # Phenotype summary
  if (nrow(gblup_data@phenotypes) > 0) {
    cat("Phenotypes:\n")
    cat("  Observations:", nrow(gblup_data@phenotypes), "\n")
    cat("  Individuals:", length(unique(gblup_data@phenotypes$GID)), "\n")
    cat("  Environments:", length(unique(gblup_data@phenotypes$ENV)), "\n")
    
    if ("trait" %in% names(gblup_data@phenotypes)) {
      trait_vals <- gblup_data@phenotypes$trait
      cat("  Trait summary:\n")
      cat("    Mean:", round(mean(trait_vals, na.rm = TRUE), 3), "\n")
      cat("    SD:", round(sd(trait_vals, na.rm = TRUE), 3), "\n")
      cat("    Range:", paste(round(range(trait_vals, na.rm = TRUE), 3), collapse = " - "), "\n")
      cat("    Missing:", sum(is.na(trait_vals)), "\n")
    }
  } else {
    cat("Phenotypes: None\n")
  }
  
  invisible(NULL)
}

#' Diagnose GBLUP Model
#'
#' Prints diagnostic information about a fitted GBLUPModel.
#'
#' @param model GBLUPModel object
#' @return Invisible NULL
#' @export
#' @examples
#' \dontrun{
#' model <- fit_gblup(gdata, K, effects = "A")
#' diagnose_model(model)
#' }
diagnose_model <- function(model) {
  
  if (!inherits(model, "GBLUPModel")) {
    stop("model must be a GBLUPModel object")
  }
  
  cat("Model Diagnostics\n")
  cat("=================\n\n")
  
  # Check convergence
  if (!is.null(model@model$convergence)) {
    cat("Convergence:", model@model$convergence, "\n")
    if (!is.null(model@model$nIter)) {
      cat("Iterations:", model@model$nIter, "\n")
    }
  }
  
  # Check variance components
  cat("\nVariance Components:\n")
  print(model@varcomp)
  
  # Check for issues
  cat("\nDiagnostic Checks:\n")
  
  # Negative variances
  if (any(model@varcomp$Variance < 0)) {
    cat("  WARNING: Negative variance components detected\n")
  } else {
    cat("  OK: All variance components non-negative\n")
  }
  
  # Zero variances
  if (any(model@varcomp$Variance == 0)) {
    cat("  WARNING: Zero variance components detected\n")
  }
  
  # GEBV statistics
  cat("\nGEBV Statistics:\n")
  gebv_vals <- model@gebv$GEBV
  cat("  N:", sum(!is.na(gebv_vals)), "\n")
  cat("  Mean:", round(mean(gebv_vals, na.rm = TRUE), 3), "\n")
  cat("  SD:", round(sd(gebv_vals, na.rm = TRUE), 3), "\n")
  cat("  Range:", paste(round(range(gebv_vals, na.rm = TRUE), 3), collapse = " to "), "\n")
  cat("  Missing:", sum(is.na(gebv_vals)), "\n")
  
  # Check GEBV variance
  if (sd(gebv_vals, na.rm = TRUE) < 1e-10) {
    cat("  WARNING: GEBV has near-zero variance\n")
  }
  
  # Observed statistics
  cat("\nObserved Statistics:\n")
  obs_vals <- model@gebv$observed
  cat("  N:", sum(!is.na(obs_vals)), "\n")
  cat("  Mean:", round(mean(obs_vals, na.rm = TRUE), 3), "\n")
  cat("  SD:", round(sd(obs_vals, na.rm = TRUE), 3), "\n")
  cat("  Range:", paste(round(range(obs_vals, na.rm = TRUE), 3), collapse = " to "), "\n")
  
  # Correlation
  complete <- complete.cases(gebv_vals, obs_vals)
  if (sum(complete) >= 2) {
    cor_val <- cor(gebv_vals[complete], obs_vals[complete])
    cat("\nPrediction Accuracy (correlation):", round(cor_val, 3), "\n")
  }
  
  invisible(NULL)
}

#' Check Data Compatibility
#'
#' Verifies that phenotype GIDs match genotype rownames and K matrix dimensions.
#'
#' @param gdata GBLUPData object
#' @param K_matrices Optional list of relationship matrices
#' @return Invisible TRUE if compatible, otherwise stops with error
#' @export
check_data_compatibility <- function(gdata, K_matrices = NULL) {
  
  if (!inherits(gdata, "GBLUPData")) {
    stop("gdata must be a GBLUPData object")
  }
  
  cat("Checking data compatibility...\n")
  
  # Get IDs from different sources
  geno_ids <- rownames(gdata@genotypes)
  pheno_ids <- unique(gdata@phenotypes$GID)
  
  # Check genotype-phenotype overlap
  overlap <- intersect(geno_ids, pheno_ids)
  
  cat("  Genotype IDs:", length(geno_ids), "\n")
  cat("  Phenotype IDs:", length(pheno_ids), "\n")
  cat("  Overlap:", length(overlap), "\n")
  
  if (length(overlap) == 0) {
    stop("No overlap between genotype and phenotype IDs!")
  }
  
  # Check K matrices if provided
  if (!is.null(K_matrices)) {
    for (nm in names(K_matrices)) {
      K <- K_matrices[[nm]]
      k_ids <- rownames(K)
      
      if (is.null(k_ids)) {
        cat("  WARNING: K matrix '", nm, "' has no rownames\n")
      } else {
        k_overlap <- length(intersect(k_ids, overlap))
        cat("  K matrix '", nm, "' overlap:", k_overlap, "\n")
        
        if (k_overlap == 0) {
          stop("K matrix '", nm, "' has no overlap with data IDs!")
        }
      }
    }
  }
  
  cat("Data compatibility: OK\n")
  invisible(TRUE)
}
#' Calculate True Prediction Accuracy from Predictive Ability
#'
#' Converts predictive ability (correlation between GEBV and phenotype) to 
#' true prediction accuracy by adjusting for heritability.
#'
#' @details
#' Predictive ability (r) is the correlation between predicted GEBVs and observed
#' phenotypes: r = cor(GEBV, y).
#'
#' True prediction accuracy (rA) is the correlation between predicted and true
#' breeding values: rA = r / sqrt(h2)
#'
#' \strong{IMPORTANT}: The h2 parameter must be the TRUE heritability, not an 
#' estimate from the training data. Using estimated h2 (especially from small
#' training sets) will produce unreliable and potentially misleading accuracy
#' values. Use this function only when:
#' \itemize{
#'   \item You know the true h2 from simulation parameters
#'   \item You have a reliable h2 estimate from a large independent dataset
#'   \item You want to make theoretical comparisons
#' }
#'
#' @param predictive_ability Correlation between GEBV and observed phenotype
#' @param h2 TRUE narrow-sense heritability (NOT estimated from training data)
#'
#' @return True prediction accuracy (correlation between predicted and true breeding values)
#' @export
#' @examples
#' \dontrun{
#' # CV returned predictive ability of 0.5
#' # TRUE heritability (from simulation) is 0.6
#' true_accuracy <- calc_prediction_accuracy(0.5, h2 = 0.6)
#' # true_accuracy = 0.5 / sqrt(0.6) = 0.645
#' 
#' # DO NOT use estimated h2 from training data - results will be unreliable
#' }
calc_prediction_accuracy <- function(predictive_ability, h2) {
  if (h2 <= 0 || h2 > 1) {
    stop("h2 must be between 0 and 1 (exclusive of 0)")
  }
  
  if (any(abs(predictive_ability) > 1, na.rm = TRUE)) {
    warning("Predictive ability should be between -1 and 1")
  }
  
  accuracy <- predictive_ability / sqrt(h2)
  
  # Cap at 1 (can exceed due to estimation error)
  accuracy[accuracy > 1] <- 1
  accuracy[accuracy < -1] <- -1
  
  return(accuracy)
}

#' Expected Prediction Accuracy (Theoretical)
#'
#' Calculates the theoretically expected prediction accuracy based on
#' training population size, heritability, and effective number of markers.
#' Uses the formula from Daetwyler et al. (2008).
#'
#' @details
#' The formula is: r = sqrt(N * h2 / (N * h2 + Me))
#'
#' Where:
#' \itemize{
#'   \item N = training population size
#'   \item h2 = narrow-sense heritability
#'   \item Me = effective number of independent chromosome segments
#' }
#'
#' Me can be approximated as: Me ≈ 2 * Ne * L where Ne is effective population
#' size and L is genome length in Morgans. For most crop species, Me ranges
#' from 100-1000. A common approximation is Me ≈ number of markers / 10 to 
#' number of markers / 100, depending on LD.
#'
#' @param n_train Training population size
#' @param h2 Narrow-sense heritability
#' @param Me Effective number of independent chromosome segments. If NULL,
#'   uses n_markers/50 as a rough approximation.
#' @param n_markers Number of markers (used to estimate Me if not provided)
#'
#' @return Expected prediction accuracy (theoretical upper bound)
#' @export
#' @references
#' Daetwyler HD, Villanueva B, Woolliams JA (2008). Accuracy of predicting the
#' genetic risk of disease using a genome-wide approach. PLoS ONE 3:e3395.
#'
#' @examples
#' # With 200 training individuals, h2=0.5, and Me=500
#' expected_accuracy(n_train = 200, h2 = 0.5, Me = 500)
#' # Returns: 0.41
#'
#' # With smaller Me (higher LD), accuracy is higher
#' expected_accuracy(n_train = 200, h2 = 0.5, Me = 100)
#' # Returns: 0.71
#'
#' # Estimate Me from markers (Me ≈ n_markers/50)
#' expected_accuracy(n_train = 200, h2 = 0.5, n_markers = 10000)
#' # Returns: 0.41 (Me = 200)
expected_accuracy <- function(n_train, h2, Me = NULL, n_markers = NULL) {
  
  if (h2 <= 0 || h2 > 1) {
    stop("h2 must be between 0 (exclusive) and 1")
  }
  
  if (n_train < 10) {
    warning("Training population size < 10 may give unreliable estimates")
  }
  
  # Estimate Me if not provided
  if (is.null(Me)) {
    if (is.null(n_markers)) {
      stop("Either Me or n_markers must be provided")
    }
    # Rough approximation: Me ≈ markers / 50 (assumes moderate LD)
    Me <- n_markers / 50
  }
  
  # Daetwyler formula
  accuracy <- sqrt(n_train * h2 / (n_train * h2 + Me))
  
  return(accuracy)
}

#' Calculate True Prediction Accuracy from Simulation Results
#'
#' For simulation studies where true breeding values (TBV) are known,
#' calculates the correlation between predicted GEBVs and true breeding values.
#'
#' @details
#' This function should only be used with data from \code{simulate_phenotypes()},
#' which stores the true breeding values in the TBV column.
#'
#' In contrast to predictive ability (correlation with phenotype), true prediction
#' accuracy measures how well we predict the actual genetic merit of individuals.
#'
#' @param model GBLUPModel object fitted to simulated data
#' @param gdata GBLUPData object with TBV column from simulate_phenotypes()
#'
#' @return Named vector with:
#'   \item{true_accuracy}{Correlation between GEBV and true additive breeding value}
#'   \item{predictive_ability}{Correlation between GEBV and phenotype}
#'   \item{theoretical_ratio}{Ratio should approximate sqrt(h2)}
#' @export
#' @examples
#' \dontrun{
#' # Simulate data
#' gdata <- simulate_genotypes(n_ind = 200, n_snp = 1000)
#' gdata <- simulate_phenotypes(gdata, h2 = 0.6)
#' K <- calc_relationship_matrices(gdata, matrices = "A")
#' model <- fit_gblup(gdata, K, effects = "A")
#' 
#' # Calculate true accuracy (only works with simulated data)
#' acc <- calc_true_accuracy(model, gdata)
#' print(acc)
#' # true_accuracy: 0.72, predictive_ability: 0.56, ratio: 0.78 (should be ~sqrt(0.6)=0.77)
#' }
calc_true_accuracy <- function(model, gdata) {
  
  if (!inherits(model, "GBLUPModel")) {
    stop("model must be a GBLUPModel object")
  }
  
  if (!inherits(gdata, "GBLUPData")) {
    stop("gdata must be a GBLUPData object")
  }
  
  # Check if TBV exists (only in simulated data)
  if (!"TBV" %in% names(gdata@phenotypes)) {
    stop("TBV column not found. This function only works with data from simulate_phenotypes().")
  }
  
  # Get GEBVs
  gebv_df <- model@gebv
  gebv_df$GID <- as.character(gebv_df$GID)
  
  # Get true breeding values (average across environments for each individual)
  pheno <- gdata@phenotypes
  pheno$GID <- as.character(pheno$GID)
  tbv_by_gid <- aggregate(cbind(TBV, trait) ~ GID, data = pheno, FUN = mean, na.rm = TRUE)
  
  # Merge
  merged <- merge(gebv_df, tbv_by_gid, by = "GID")
  
  if (nrow(merged) == 0) {
    stop("No matching GIDs between model GEBVs and phenotype data")
  }
  
  # Calculate correlations
  valid <- complete.cases(merged[, c("GEBV", "TBV", "trait")])
  
  true_accuracy <- cor(merged$GEBV[valid], merged$TBV[valid])
  predictive_ability <- cor(merged$GEBV[valid], merged$trait[valid])
  
  # The ratio should approximate sqrt(h2)
  ratio <- predictive_ability / true_accuracy
  
  return(c(
    true_accuracy = true_accuracy,
    predictive_ability = predictive_ability,
    theoretical_ratio = ratio,
    n = sum(valid)
  ))
}
