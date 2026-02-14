#' GBLUPData Class
#'
#' S4 class to store genotype and phenotype data for GBLUP analysis
#'
#' @slot genotypes Matrix of genotype data (individuals x markers), coded as -1, 0, 1
#' @slot phenotypes Data frame with columns: GID (genotype ID), ENV (environment), trait (phenotypic value)
#' @slot maf Numeric vector of minor allele frequencies
#' @slot metadata List of additional metadata
#'
#' @export
setClass("GBLUPData",
         slots = c(
           genotypes = "matrix",
           phenotypes = "data.frame",
           maf = "numeric",
           metadata = "list"
         ),
         prototype = list(
           metadata = list()
         ),
         validity = function(object) {
           errors <- character()
           
           # Check phenotypes has required columns
           required_cols <- c("GID", "ENV", "trait")
           if (!all(required_cols %in% names(object@phenotypes))) {
             errors <- c(errors, "phenotypes must have columns: GID, ENV, trait")
           }
           
           # Check maf length matches number of SNPs
           if (length(object@maf) != ncol(object@genotypes)) {
             errors <- c(errors, "maf length must match number of SNPs")
           }
           
           # Robust genotype coding check
           allowed_genotype_values <- c(-1, 0, 1)
           unique_geno_values <- unique(as.vector(object@genotypes))
           if (!all(unique_geno_values %in% allowed_genotype_values)) {
             not_allowed <- setdiff(unique_geno_values, allowed_genotype_values)
             errors <- c(errors, paste0(
               "genotypes must be coded as exactly -1, 0, 1; found invalid values: ",
               paste(sort(not_allowed), collapse = ", ")
             ))
           }
           
           if (length(errors) == 0) TRUE else errors
         }
)

#' CVResult Class
#'
#' S4 class to store cross-validation results
#'
#' @slot predictions Data frame with all predictions across folds
#' @slot metrics Data frame with metrics by fold
#' @slot fold_results List of detailed results by fold
#' @slot scheme CV scheme used
#' @slot n_folds Number of folds
#' @slot n_reps Number of replications
#' @slot effects Effects included in model
#'
#' @export
setClass("CVResult",
         slots = c(
           predictions = "data.frame",
           metrics = "data.frame",
           fold_results = "list",
           scheme = "character",
           n_folds = "numeric",
           n_reps = "numeric",
           effects = "character"
         ),
         prototype = list(
           scheme = "CV1",
           n_folds = 5,
           n_reps = 1,
           effects = "A"
         )
)

#' GWASResult Class
#'
#' S4 class to store GWAS results
#'
#' @slot results Data frame with marker test results
#' @slot lambda Genomic inflation factor
#' @slot n_markers Number of markers tested
#' @slot min_MAF Minimum MAF used
#' @slot method Test method used
#'
#' @export
setClass("GWASResult",
         slots = c(
           results = "data.frame",
           lambda = "numeric",
           n_markers = "numeric",
           min_MAF = "numeric",
           method = "character"
         ),
         prototype = list(
           lambda = NA_real_,
           n_markers = 0,
           min_MAF = 0.05,
           method = "Wald"
         )
)

#' GBLUPModel Class
#'
#' S4 class to store results from GBLUP model fitting
#'
#' @slot model Fitted model object from internal MME solver
#' @slot data Original GBLUPData object
#' @slot gebv Data frame with genomic estimated breeding values
#' @slot varcomp Data frame with variance components (columns: Component, Variance)
#' @slot description Character string describing the model
#' @slot effects Character vector of effects included in the model
#'
#' @export
setClass("GBLUPModel",
         slots = c(
           model = "ANY",
           data = "GBLUPData",
           gebv = "data.frame",
           varcomp = "data.frame",
           description = "character",
           effects = "character"
         ),
         prototype = list(
           description = "GBLUP Model",
           effects = "A"
         )
)

#' Show Methods for CompGEBLUP Objects
#'
#' Display summary information for S4 classes
#'
#' @param object An S4 object (GBLUPData, GBLUPModel, CVResult, or GWASResult)
#' @name show-methods
#' @rdname show-methods
NULL

#' @rdname show-methods
#' @export
setMethod("show", "GBLUPData", function(object) {
  cat("GBLUPData object\n")
  cat("================\n")
  cat("Genotypes:", nrow(object@genotypes), "individuals x", 
      ncol(object@genotypes), "SNPs\n")
  cat("Phenotypes:", nrow(object@phenotypes), "observations\n")
  cat("Environments:", length(unique(object@phenotypes$ENV)), "\n")
  cat("MAF range:", round(range(object@maf), 3), "\n")
  cat("Missing data:", 
      round(100 * sum(is.na(object@genotypes)) / length(object@genotypes), 2), "%\n")
})

#' @rdname show-methods
#' @export
setMethod("show", "GBLUPModel", function(object) {
  cat("GBLUPModel object\n")
  cat("=================\n")
  
  if (nchar(object@description) > 0 && object@description != "GBLUP Model") {
    cat("Model:", object@description, "\n")
  } else {
    cat("Model type: GBLUP\n")
  }
  
  cat("Effects:", paste(object@effects, collapse = " + "), "\n")
  cat("Individuals:", length(unique(object@gebv$GID)), "\n")
  if ("ENV" %in% names(object@gebv)) {
    cat("Environments:", length(unique(object@gebv$ENV)), "\n")
  }
  
  cat("\nVariance components:\n")
  print(object@varcomp, row.names = FALSE)
  
  # Calculate and show heritability if possible
  if ("A" %in% object@effects && nrow(object@varcomp) > 0) {
    va <- object@varcomp$Variance[grep("^A$|^A\\.|GID", object@varcomp$Component)]
    vr <- object@varcomp$Variance[grep("Residual", object@varcomp$Component, ignore.case = TRUE)]
    if (length(va) > 0 && length(vr) > 0) {
      total_var <- sum(object@varcomp$Variance)
      h2 <- sum(va) / total_var
      cat("\nNarrow-sense heritability (h2):", round(h2, 3), "\n")
    }
  }
})

#' @rdname show-methods
#' @export
setMethod("show", "CVResult", function(object) {
  cat("Cross-Validation Result\n")
  cat("=======================\n")
  cat("Scheme:", object@scheme, "\n")
  cat("Number of folds:", object@n_folds, "\n")
  cat("Number of replications:", object@n_reps, "\n")
  cat("Effects included:", paste(object@effects, collapse = ", "), "\n")
  cat("\n")
  
  if (nrow(object@metrics) > 0) {
    cat("Summary Metrics:\n")
    
    # Predictive ability - THE primary metric
    if ("predictive_ability" %in% names(object@metrics)) {
      pa_vals <- object@metrics$predictive_ability[!is.na(object@metrics$predictive_ability)]
      if (length(pa_vals) > 0) {
        cat("  Predictive ability:  ", round(mean(pa_vals), 4), 
            " (SD:", round(sd(pa_vals), 4), ")\n", sep = "")
      }
    }
    
    # MSE
    if ("mse" %in% names(object@metrics)) {
      cat("  Mean MSE:            ", round(mean(object@metrics$mse, na.rm = TRUE), 4), "\n", sep = "")
    }
    
    cat("  Total predictions:   ", nrow(object@predictions), "\n", sep = "")
    
    cat("\nNote: Predictive ability = cor(GEBV, phenotype)\n")
    cat("      To convert to accuracy: calc_prediction_accuracy(r, h2_true)\n")
    cat("      For theoretical max: expected_accuracy(n_train, h2, Me)\n")
  } else {
    cat("No valid results available\n")
  }
})

#' @rdname show-methods
#' @export
setMethod("show", "GWASResult", function(object) {
  cat("GWAS Result\n")
  cat("===========\n")
  cat("Number of markers tested:", object@n_markers, "\n")
  cat("Minimum MAF:", object@min_MAF, "\n")
  cat("Method:", object@method, "\n")
  cat("Genomic inflation (lambda):", round(object@lambda, 3), "\n")
  cat("\n")
  
  # Summary of results
  bonf_thresh <- 0.05 / object@n_markers
  n_sig <- sum(object@results$p_value < bonf_thresh, na.rm = TRUE)
  n_sug <- sum(object@results$p_value < 1e-5, na.rm = TRUE)
  
  cat("Significant markers (Bonferroni):", n_sig, "\n")
  cat("Suggestive markers (p < 1e-5):", n_sug, "\n")
  cat("\n")
  
  # Top markers
  cat("Top 10 markers:\n")
  top <- get_top_markers(object, n = 10)
  print(top[, c("marker", "beta", "p_value", "log10p")])
})

#' Calculate Predictive Ability (Correlation)
#'
#' Calculates correlation between observed phenotypes and predicted GEBVs (predictive ability).
#'
#' @param object A GBLUPModel object
#' @return Numeric correlation coefficient between observed and predicted values
#' @rdname calculate_accuracy
#' @export
setGeneric("calculate_accuracy", function(object) {
  standardGeneric("calculate_accuracy")
})

#' @rdname calculate_accuracy
#' @export
setMethod("calculate_accuracy", "GBLUPModel", function(object) {
  gebv <- object@gebv
  complete_cases <- complete.cases(gebv[, c("GEBV", "observed")])
  
  if (sum(complete_cases) == 0) {
    warning("No complete cases for accuracy calculation")
    return(NA)
  }
  
  cor(gebv$GEBV[complete_cases], 
      gebv$observed[complete_cases])
})

#' Summary method for GBLUPModel
#' @param object A GBLUPModel object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the object
#' @export
setMethod("summary", "GBLUPModel", function(object, ...) {
  cat("GBLUP Model Summary\n")
  cat("===================\n\n")
  
  show(object)
  
  # GEBV Statistics
  cat("\nGEBV Statistics:\n")
  cat("  Mean:", round(mean(object@gebv$GEBV, na.rm = TRUE), 3), "\n")
  cat("  SD:", round(sd(object@gebv$GEBV, na.rm = TRUE), 3), "\n")
  cat("  Range:", paste(round(range(object@gebv$GEBV, na.rm = TRUE), 3), collapse = " to "), "\n\n")
  
  # Prediction Accuracy
  acc <- calculate_accuracy(object)
  cat("Prediction Accuracy (correlation):", round(acc, 3), "\n")
  
  invisible(object)
})

#' Summary Method for CVResult
#'
#' @param object CVResult object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the object
#' @export
setMethod("summary", "CVResult", function(object, ...) {
  cat("Cross-Validation Summary\n")
  cat("========================\n\n")
  
  show(object)
  
  if (nrow(object@metrics) > 0) {
    cat("\nMetrics by Fold:\n")
    print(object@metrics)
    
    if (object@n_reps > 1) {
      cat("\nMetrics by Replication:\n")
      rep_summary <- aggregate(cbind(predictive_ability, mse) ~ replication, 
                               data = object@metrics,
                               FUN = function(x) c(mean = mean(x), sd = sd(x)))
      print(rep_summary)
    }
  }
  
  invisible(object)
})

#' Summary Method for GWASResult
#'
#' @param object GWASResult object
#' @param ... Additional arguments (ignored)
#' @return Invisibly returns the object
#' @export
setMethod("summary", "GWASResult", function(object, ...) {
  show(object)
  
  cat("\nDistribution of p-values:\n")
  print(summary(object@results$p_value))
  
  cat("\nDistribution of effect sizes:\n")
  print(summary(object@results$beta))
  
  invisible(object)
})

#' Predict method for GBLUPModel
#' @param object A GBLUPModel object
#' @param ... Additional arguments (not used)
#' @return Data frame with predictions
#' @export
setMethod("predict", "GBLUPModel", function(object, ...) {
  return(object@gebv[, c("GID", "ENV", "GEBV")])
})

#' Extract GEBV from Model
#'
#' @param model GBLUPModel object
#' @return Data frame with GID, ENV, GEBV, and observed values
#' @export
#' @examples
#' \dontrun{
#' model <- fit_gblup(gdata, K, effects = "A")
#' gebv <- extract_gebv(model)
#' head(gebv)
#' }
extract_gebv <- function(model) {
  if (!inherits(model, "GBLUPModel")) {
    stop("model must be a GBLUPModel object")
  }
  return(model@gebv)
}

#' Extract Variance Components
#'
#' @param model GBLUPModel object
#' @return Data frame with variance components (columns: Component, Variance)
#' @export
#' @examples
#' \dontrun{
#' model <- fit_gblup(gdata, K, effects = "A")
#' vc <- extract_varcomp(model)
#' print(vc)
#' }
extract_varcomp <- function(model) {
  if (!inherits(model, "GBLUPModel")) {
    stop("model must be a GBLUPModel object")
  }
  return(model@varcomp)
}