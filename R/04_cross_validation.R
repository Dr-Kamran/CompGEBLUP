#' Cross-Validation for GBLUP Models
#'
#' Performs k-fold cross-validation for genomic prediction models.
#'
#' @details
#' The returned metrics include:
#' \itemize{
#'   \item \strong{predictive_ability}: Correlation between predicted GEBVs and observed 
#'     phenotypes, cor(GEBV, y). This is THE standard metric reported in genomic prediction
#'     literature. A value of ~0.4 is typical for small training populations.
#'   \item \strong{MSE}: Mean squared error of predictions
#' }
#'
#' Note: "Predictive ability" (correlation with phenotype) differs from "prediction accuracy"
#' (correlation with true breeding value). The theoretical relationship is:
#' accuracy = predictive_ability / sqrt(h²). However, this conversion requires knowing the
#' TRUE heritability. Using estimated h² from training data produces unreliable results.
#' 
#' For theoretical comparisons, use \code{expected_accuracy()} to estimate the theoretical
#' maximum accuracy based on training population size and known h².
#'
#' Cross-validation schemes:
#' \itemize{
#'   \item \strong{CV1}: Random k-fold - standard evaluation of genomic prediction
#'   \item \strong{CV2}: Leave-one-environment-out - evaluates prediction to new environments
#'   \item \strong{CV0}: Within-environment - evaluates prediction within each environment
#' }
#'
#' @param gdata GBLUPData object
#' @param K_matrices List of relationship matrices (used only when recompute_K = FALSE)
#' @param effects Character vector of random effects to include
#' @param scheme Cross-validation scheme: "CV1" (random), "CV2" (leave-one-environment-out),
#'   or "CV0" (within-environment)
#' @param n_folds Number of folds for CV0 and CV1 schemes
#' @param n_reps Number of replications
#' @param ridge Ridge parameter for numerical stability (scalar, applied to all matrices).
#'   Overridden component-by-component by ridge_per_matrix when that argument is provided.
#' @param ridge_per_matrix Optional named numeric vector of per-matrix ridge values.
#'   Names must match matrix types: "A", "D", "AA", "AE", "DE", "AAE".
#'   Components not listed fall back to the scalar \code{ridge}.
#'   Example: \code{ridge_per_matrix = c(A = 0.001, D = 0.05, AA = 0.05)}.
#'   Only used when \code{recompute_K = TRUE}. This is the correct way to apply
#'   differential regularisation (e.g. stronger ridge on ill-conditioned D or AA)
#'   without weakening the baseline A matrix regularisation.
#' @param recompute_K Logical, if TRUE (default) recomputes K matrix for each fold using
#'   only training individuals. This prevents data leakage from test genotypes but is slower.
#'   If FALSE, uses pre-computed K_matrices for all folds (faster but allows minor leakage).
#' @param verbose Logical, print progress messages
#' @param fixed Optional fixed effects specification. Can be NULL (intercept only, default),
#'   a one-sided formula (e.g., ~ TESTER), or a pre-built design matrix. Fixed effects
#'   are applied to the training data in each fold.
#'
#' @return CVResult object containing:
#'   \item{metrics}{Data frame with predictive_ability and MSE per fold}
#'   \item{predictions}{Data frame with predicted vs observed values}
#' @export
#' @examples
#' \dontrun{
#' gdata <- simulate_genotypes(n_ind = 100, n_snp = 500)
#' gdata <- simulate_phenotypes(gdata, n_env = 3, h2 = 0.6)
#' K <- calc_relationship_matrices(gdata, matrices = "A")
#' 
#' # Run CV - reports predictive ability (the standard GP metric)
#' cv_res <- cv_gblup(gdata, K, effects = "A", scheme = "CV1", n_folds = 5)
#' 
#' # With fixed effects
#' cv_res_fixed <- cv_gblup(gdata, K, effects = "A", fixed = ~ TESTER, 
#'                          scheme = "CV1", n_folds = 5)
#' 
#' # Compare to theoretical maximum (requires knowing true h² and Me)
#' # With Me=500, theoretical max is ~0.30; with Me=100, it's ~0.57
#' expected_accuracy(n_train = 80, h2 = 0.6, Me = 100)
#' 
#' # Convert to accuracy only if you know the TRUE h² (e.g., from simulation)
#' calc_prediction_accuracy(cv_res@metrics$predictive_ability, h2 = 0.6)
#' }
cv_gblup <- function(gdata, K_matrices = NULL, effects = "A",
                     scheme = c("CV1", "CV2", "CV0"), 
                     n_folds = 5, n_reps = 1,
                     ridge = 0.001, ridge_per_matrix = NULL,
                     recompute_K = TRUE, verbose = TRUE,
                     fixed = NULL) {
  
  # Match arguments
  scheme <- match.arg(scheme)
  
  if (!inherits(gdata, "GBLUPData")) {
    stop("gdata must be a GBLUPData object")
  }
  
  # Validate K_matrices when not recomputing
  if (!recompute_K && is.null(K_matrices)) {
    stop("K_matrices must be provided when recompute_K = FALSE")
  }
  
  # Warn about data leakage when using pre-computed K
  if (!recompute_K) {
    warning("recompute_K = FALSE: Using pre-computed K matrices may cause ",
            "data leakage (test genotypes influence K). ",
            "Prediction accuracies may be inflated by 5-15%. ",
            "Set recompute_K = TRUE for unbiased CV (slower but correct).",
            call. = FALSE)
  }
  
  stopifnot("n_folds must be at least 2" = n_folds >= 2)
  stopifnot("n_reps must be at least 1" = n_reps >= 1)
  
  # Validate and normalise ridge_per_matrix
  # BUG FIX (v2.2.1): allow per-matrix ridge values so D and AA can be
  # regularised more strongly than A, matching the full model behaviour.
  valid_matrix_types <- c("A", "D", "AA", "AE", "DE", "AAE")
  if (!is.null(ridge_per_matrix)) {
    if (!is.numeric(ridge_per_matrix) || is.null(names(ridge_per_matrix))) {
      stop("ridge_per_matrix must be a named numeric vector, e.g. c(A=0.001, D=0.05, AA=0.05)")
    }
    bad_names <- setdiff(names(ridge_per_matrix), valid_matrix_types)
    if (length(bad_names) > 0) {
      stop("ridge_per_matrix has unknown matrix type(s): ",
           paste(bad_names, collapse = ", "),
           ". Valid names: ", paste(valid_matrix_types, collapse = ", "))
    }
    if (any(ridge_per_matrix < 0)) {
      stop("All values in ridge_per_matrix must be >= 0")
    }
    if (verbose && recompute_K) {
      cat("Differential ridge (ridge_per_matrix):\n")
      for (nm in names(ridge_per_matrix)) {
        cat(sprintf("  %s: %.4f\n", nm, ridge_per_matrix[nm]))
      }
      unlisted <- setdiff(effects[effects %in% valid_matrix_types], names(ridge_per_matrix))
      if (length(unlisted) > 0) {
        cat(sprintf("  %s: %.4f (default ridge)\n",
                    paste(unlisted, collapse = ", "), ridge))
      }
    }
  }
  
  # Get phenotype data
  pheno <- gdata@phenotypes
  
  # Ensure character types
  pheno$GID <- as.character(pheno$GID)
  pheno$ENV <- as.character(pheno$ENV)
  
  # Remove any existing NA values
  pheno_complete <- pheno[complete.cases(pheno$trait), ]
  
  if (nrow(pheno_complete) < n_folds) {
    stop("Not enough observations for ", n_folds, " folds. Have ", 
         nrow(pheno_complete), " complete observations.")
  }
  
  # Update gdata with complete phenotypes only
  gdata@phenotypes <- pheno_complete
  
  if (verbose && recompute_K) {
    cat("Using strict CV: K matrix recomputed per fold (no data leakage)\n")
  } else if (verbose) {
    cat("Using fast CV: K matrix computed once (minor genotype leakage)\n")
  }
  
  # Initialize results storage
  all_results <- list()
  fold_errors <- list()
  
  # Track timing for progress estimates
  start_time <- Sys.time()
  total_folds_all_reps <- 0
  completed_folds <- 0
  
  # Perform replications
  for (rep in 1:n_reps) {
    
    if (verbose && n_reps > 1) {
      cat("Replication", rep, "of", n_reps, "\n")
    }
    
    # Create folds based on scheme
    folds <- switch(scheme,
                    "CV1" = create_random_folds(pheno_complete, n_folds),
                    "CV2" = create_env_folds(pheno_complete),
                    "CV0" = create_within_env_folds(pheno_complete, n_folds)
    )
    
    # Count total folds on first rep
    if (rep == 1) {
      total_folds_all_reps <- length(folds) * n_reps
    }
    
    # Perform cross-validation
    rep_results <- lapply(seq_along(folds), function(i) {
      
      if (verbose) {
        completed_folds <<- completed_folds + 1
        elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
        if (completed_folds > 1) {
          avg_time_per_fold <- elapsed / (completed_folds - 1)
          remaining_folds <- total_folds_all_reps - completed_folds
          eta_secs <- avg_time_per_fold * remaining_folds
          eta_str <- if (eta_secs > 60) {
            sprintf("%.1f min", eta_secs / 60)
          } else {
            sprintf("%.0f sec", eta_secs)
          }
          cat("  Fold", i, "of", length(folds), "- ETA:", eta_str, "\n")
        } else {
          cat("  Fold", i, "of", length(folds), "\n")
        }
      }
      
      fold <- folds[[i]]
      
      # Fit model and predict
      result <- tryCatch({
        cv_fold(gdata, K_matrices, effects, fold, ridge, 
                ridge_per_matrix = ridge_per_matrix,
                recompute_K = recompute_K, verbose = FALSE, fixed = fixed)
      }, error = function(e) {
        if (verbose) {
          cat("    ERROR in fold", i, ":", e$message, "\n")
        }
        fold_errors[[length(fold_errors) + 1]] <<- list(
          fold = i,
          rep = rep,
          error = e$message
        )
        list(
          fold = i,
          observed = numeric(0),
          predicted = numeric(0),
          predictive_ability = NA,
          mse = NA,
          error = e$message
        )
      })
      
      result$fold <- i
      result$replication <- rep
      
      return(result)
    })
    
    all_results[[rep]] <- rep_results
  }
  
  # Report total time
  if (verbose) {
    total_time <- difftime(Sys.time(), start_time, units = "secs")
    cat("Cross-validation completed in", 
        sprintf("%.1f", as.numeric(total_time)), "seconds\n")
  }
  
  # Check if all folds failed
  all_failed <- all(sapply(unlist(all_results, recursive = FALSE), function(x) {
    is.null(x$observed) || length(x$observed) == 0
  }))
  
  if (all_failed) {
    if (length(fold_errors) > 0) {
      cat("\nFold errors:\n")
      for (err in fold_errors) {
        cat("  Rep", err$rep, "Fold", err$fold, ":", err$error, "\n")
      }
    }
    stop("All folds failed. Check your data and parameters.")
  }
  
  # Summarize results
  summary_stats <- summarize_cv_results(all_results, scheme)
  
  # Create CVResult object
  cv_result <- new("CVResult",
                   predictions = compile_predictions(all_results),
                   metrics = summary_stats$metrics,
                   fold_results = all_results,
                   scheme = scheme,
                   n_folds = if (scheme == "CV2") length(unique(pheno_complete$ENV)) else n_folds,
                   n_reps = n_reps,
                   effects = effects)
  
  if (verbose) {
    cat("\nCross-validation completed\n")
    cat("Mean predictive ability:", round(summary_stats$mean_predictive_ability, 3), "\n")
    cat("SD predictive ability:  ", round(summary_stats$sd_predictive_ability, 3), "\n")
    cat("Mean MSE:", round(summary_stats$mean_mse, 3), "\n")
  }
  
  return(cv_result)
}

#' Create Random Folds for CV1
#'
#' CV1 scheme: Random cross-validation across all environments.
#' All observations of each individual are kept together (all in train OR all in test).
#'
#' @param pheno Phenotype data frame
#' @param n_folds Number of folds
#'
#' @return List of fold indices
#' @keywords internal
create_random_folds <- function(pheno, n_folds) {
  
  # Ensure character GID
  pheno$GID <- as.character(pheno$GID)
  
  # Get unique individuals
  individuals <- unique(pheno$GID)
  n_ind <- length(individuals)
  
  if (n_ind < n_folds) {
    stop("Not enough unique individuals (", n_ind, ") for ", n_folds, " folds")
  }
  
  # Shuffle individuals
  shuffled <- sample(individuals)
  
  # Create fold assignments
  fold_assignments <- rep(1:n_folds, length.out = n_ind)
  
  folds <- lapply(1:n_folds, function(i) {
    test_ids <- shuffled[fold_assignments == i]
    list(
      # Test includes ALL observations (all environments) for test individuals
      test = which(pheno$GID %in% test_ids),
      # Train includes ALL observations (all environments) for train individuals
      train = which(!pheno$GID %in% test_ids),
      test_individuals = test_ids
    )
  })
  
  return(folds)
}

#' Create Environment Folds for CV2
#'
#' CV2 scheme: Leave-one-environment-out cross-validation.
#' Tests ability to predict in completely new environments.
#'
#' @param pheno Phenotype data frame
#'
#' @return List of fold indices
#' @keywords internal
create_env_folds <- function(pheno) {
  
  # Ensure character ENV
  pheno$ENV <- as.character(pheno$ENV)
  
  # Get unique environments
  environments <- unique(pheno$ENV)
  
  if (length(environments) < 2) {
    stop("Need at least 2 environments for CV2 scheme")
  }
  
  # Create one fold per environment
  folds <- lapply(environments, function(env) {
    list(
      test = which(pheno$ENV == env),
      train = which(pheno$ENV != env),
      env = env,
      test_individuals = unique(pheno$GID[pheno$ENV == env])
    )
  })
  
  return(folds)
}

#' Create Within-Environment Folds for CV0
#'
#' CV0 scheme: Random cross-validation within each environment.
#' 
#' CRITICAL: To prevent data leakage, training set includes ONLY
#' other individuals from the SAME environment. Test individuals
#' are completely excluded from training, even if they appear in
#' other environments.
#'
#' @param pheno Phenotype data frame
#' @param n_folds Number of folds per environment
#'
#' @return List of fold indices
#' @keywords internal
create_within_env_folds <- function(pheno, n_folds) {
  
  # Ensure character types
  pheno$GID <- as.character(pheno$GID)
  pheno$ENV <- as.character(pheno$ENV)
  
  # Get unique environments
  environments <- unique(pheno$ENV)
  
  # Pre-check: all environments must have minimum individuals
  # Require at least 10 for reliable CV estimates (5 is too few for most use cases)
  min_required <- max(10, n_folds)
  small_envs <- character()
  
  for (env in environments) {
    env_ids <- unique(pheno$GID[pheno$ENV == env])
    n_ind <- length(env_ids)
    
    if (n_ind < 10) {
      stop("Environment '", env, "' has only ", n_ind, " individuals. ",
           "CV0 requires at least 10 individuals per environment for reliable estimates. ",
           "Consider using CV1 (random across environments) or CV2 (leave-one-environment-out) instead.")
    }
    
    if (n_ind < n_folds) {
      small_envs <- c(small_envs, paste0(env, " (n=", n_ind, ")"))
    }
  }
  
  # Warn about environments that will use fewer folds
  if (length(small_envs) > 0) {
    warning("The following environments have fewer individuals than n_folds=", n_folds, 
            " and will use fewer folds: ", paste(small_envs, collapse = ", "), 
            ". Results may be less reliable.", call. = FALSE)
  }
  
  # Create folds within each environment
  all_folds <- list()
  fold_counter <- 1
  
  for (env in environments) {
    env_indices <- which(pheno$ENV == env)
    env_ids <- unique(pheno$GID[env_indices])
    n_ind <- length(env_ids)
    
    if (n_ind < n_folds) {
      n_folds_env <- n_ind  # Use leave-one-out if very small
    } else {
      n_folds_env <- n_folds
    }
    
    # Shuffle individuals within environment
    shuffled <- sample(env_ids)
    
    # Create fold assignments
    fold_assignments <- rep(1:n_folds_env, length.out = n_ind)
    
    for (i in 1:n_folds_env) {
      test_ids <- shuffled[fold_assignments == i]
      train_ids <- shuffled[fold_assignments != i]
      
      all_folds[[fold_counter]] <- list(
        # Test: only observations in this environment for test individuals
        test = which(pheno$GID %in% test_ids & pheno$ENV == env),
        # Train: only observations in this environment for OTHER individuals
        # This prevents data leakage from other environments
        train = which(pheno$GID %in% train_ids & pheno$ENV == env),
        env = env,
        test_env = env,
        test_individuals = test_ids,
        train_individuals = train_ids
      )
      fold_counter <- fold_counter + 1
    }
  }
  
  return(all_folds)
}

#' Cross-Validation Fold
#'
#' Fits model on training data and predicts test data for one fold.
#'
#' @param gdata GBLUPData object
#' @param K_matrices List of relationship matrices (pre-computed, used when recompute_K = FALSE)
#' @param effects Character vector of effects
#' @param fold List with train and test indices
#' @param ridge Ridge parameter (scalar fallback for all matrices)
#' @param ridge_per_matrix Optional named numeric vector of per-matrix ridge values
#'   (names: "A", "D", "AA", "AE", "DE", "AAE"). Overrides scalar \code{ridge}
#'   for named components. Only used when \code{recompute_K = TRUE}.
#' @param recompute_K Logical, if TRUE recomputes K from training genotypes only
#' @param verbose Logical
#' @param fixed Optional fixed effects specification (formula or matrix)
#'
#' @return List with fold results
#' @keywords internal
cv_fold <- function(gdata, K_matrices, effects, fold, ridge,
                    ridge_per_matrix = NULL,
                    recompute_K = TRUE, verbose = FALSE, fixed = NULL) {
  
  # Validate fold
  if (length(fold$train) == 0) {
    stop("Training set is empty")
  }
  
  if (length(fold$test) == 0) {
    stop("Test set is empty")
  }
  
  # Get full phenotype data
  pheno_full <- gdata@phenotypes
  pheno_full$GID <- as.character(pheno_full$GID)
  pheno_full$ENV <- as.character(pheno_full$ENV)
  
  # Get environment info
  env_names <- unique(pheno_full$ENV)
  n_env <- length(env_names)
  
  # Get genotype matrix
  M <- gdata@genotypes
  all_gids <- rownames(M)
  
  # Identify training and test individuals
  train_gids <- unique(pheno_full$GID[fold$train])
  test_gids <- unique(pheno_full$GID[fold$test])
  
  # Ensure they exist in genotype matrix
  train_gids <- train_gids[train_gids %in% all_gids]
  test_gids <- test_gids[test_gids %in% all_gids]
  
  if (length(train_gids) < 5) {
    stop("Too few training individuals with genotypes (", length(train_gids), ")")
  }
  
  # Recompute K matrices if requested (strict CV - no data leakage)
  if (recompute_K) {
    # Get training genotypes only
    M_train <- M[train_gids, , drop = FALSE]
    
    # Calculate MAF from training only
    p_train <- colMeans(M_train + 1, na.rm = TRUE) / 2
    maf_train <- pmin(p_train, 1 - p_train)
    
    # Filter markers with low MAF in training (silently)
    keep_markers <- maf_train >= 0.01
    if (sum(keep_markers) < 100) {
      # Silently use all markers if too few pass filter
      keep_markers <- rep(TRUE, ncol(M_train))
    }
    M_train_filtered <- M_train[, keep_markers, drop = FALSE]
    
    # BUG FIX (v2.2.1): resolve ridge for each matrix type individually.
    # If ridge_per_matrix is provided, named components use their own value;
    # everything else falls back to the scalar `ridge`.
    .get_ridge <- function(mat_type) {
      if (!is.null(ridge_per_matrix) && mat_type %in% names(ridge_per_matrix)) {
        ridge_per_matrix[[mat_type]]
      } else {
        ridge
      }
    }
    
    # Initialize K_matrices_fold
    K_matrices_fold <- list()
    
    # Compute A matrix if needed
    if (any(c("A", "AE", "AA", "AAE") %in% effects)) {
      K_train <- build_A_matrix(M_train_filtered, verbose = FALSE)
      # BUG FIX (v2.3.1): Use ridge= (diagonal) not tol= (eigenvalue floor).
      K_train <- make_positive_definite(K_train, ridge = .get_ridge("A"))
      K_matrices_fold$A <- K_train
    }
    
    # Compute D matrix if needed
    if (any(c("D", "DE") %in% effects)) {
      D_train <- build_D_matrix(M_train_filtered, verbose = FALSE)
      D_train <- make_positive_definite(D_train, ridge = .get_ridge("D"))
      K_matrices_fold$D <- D_train
    }
    
    # Compute AA (epistasis) matrix if needed
    if (any(c("AA", "AAE") %in% effects)) {
      if (!"A" %in% names(K_matrices_fold)) {
        K_train <- build_A_matrix(M_train_filtered)
        K_train <- make_positive_definite(K_train, ridge = .get_ridge("A"))
      } else {
        K_train <- K_matrices_fold$A
      }
      # BUG FIX (v2.3.1): Use build_E_matrix for consistency with full-data build.
      # Raw Hadamard K*K skips any normalisation that build_E_matrix applies.
      AA_train <- tryCatch(
        build_E_matrix(M = M_train_filtered, type = "A#A", A = K_train,
                       min.MAF = 0, use.cpp = TRUE),
        error = function(e) K_train * K_train  # fallback to Hadamard
      )
      AA_train <- make_positive_definite(AA_train, ridge = .get_ridge("AA"))
      K_matrices_fold$AA <- AA_train
    }
    
    # Compute GxE matrices if needed
    if ("AE" %in% effects) {
      A_base <- K_matrices_fold$A
      K_matrices_fold$AE <- build_GE_matrix(A_base, n_env = n_env, env_names = env_names)
    }
    
    if ("DE" %in% effects) {
      D_base <- K_matrices_fold$D
      K_matrices_fold$DE <- build_GE_matrix(D_base, n_env = n_env, env_names = env_names)
    }
    
    if ("AAE" %in% effects) {
      AA_base <- K_matrices_fold$AA
      K_matrices_fold$AAE <- build_GE_matrix(AA_base, n_env = n_env, env_names = env_names)
    }
    
    # Pass through custom-effect K matrices that aren't rebuilt per fold.
    # These are user-provided kernels (e.g. "RN", "K_comb") that don't have
    # a built-in rebuild path. When recompute_K=TRUE, only A/D/AA/AE/DE/AAE
    # are recomputed from genotypes; custom kernels are taken from the
    # original K_matrices (which may be slightly stale for the fold's subset,
    # but this is the best we can do without user-provided rebuild logic).
    builtin_k_types <- c("A", "D", "AA", "AE", "DE", "AAE")
    for (eff_name in names(K_matrices)) {
      if (!(eff_name %in% builtin_k_types) && !(eff_name %in% names(K_matrices_fold))) {
        K_matrices_fold[[eff_name]] <- K_matrices[[eff_name]]
      }
    }
    
  } else {
    # Use pre-computed K matrices
    K_matrices_fold <- K_matrices
  }
  
  # Create masked phenotypes (set test to NA)
  pheno_masked <- pheno_full
  pheno_masked$trait[fold$test] <- NA
  
  # Check that training set has data
  n_train_obs <- sum(!is.na(pheno_masked$trait))
  if (n_train_obs < 5) {
    stop("Too few training observations (", n_train_obs, ")")
  }
  
  # Create training GBLUPData object
  if (recompute_K) {
    # Use only training individuals
    M_for_train <- M[train_gids, , drop = FALSE]
    pheno_train_only <- pheno_masked[pheno_masked$GID %in% train_gids, , drop = FALSE]
    
    # Recalculate MAF for subset
    p_subset <- colMeans(M_for_train + 1, na.rm = TRUE) / 2
    maf_subset <- pmin(p_subset, 1 - p_subset)
    
    gdata_train <- new("GBLUPData",
                       phenotypes = pheno_train_only,
                       genotypes = M_for_train,
                       maf = maf_subset,
                       metadata = list())
  } else {
    gdata_train <- new("GBLUPData",
                       phenotypes = pheno_masked,
                       genotypes = gdata@genotypes,
                       maf = gdata@maf,
                       metadata = list())
  }
  
  # Fit model on training data
  # K matrices are already regularized by make_positive_definite or taken from pre-computed
  model <- fit_gblup(gdata_train, K_matrices_fold, effects, 
                     ridge = ridge, verbose = FALSE,
                     K_already_regularized = TRUE, fixed = fixed)
  
  # Get GEBVs from training
  gebv_train <- model@gebv
  gebv_train$GID <- as.character(gebv_train$GID)
  
  # Predict test individuals
  if (recompute_K) {
    # Use proper BLUP equation: u_test = K[test,train] * K[train,train]^-1 * u_train
    # First, compute K between test and training from genotypes
    M_test <- M[test_gids, , drop = FALSE]
    M_train <- M[train_gids, , drop = FALSE]
    
    # Recalculate MAF and filter markers (same as during K computation)
    p_train_pred <- colMeans(M_train + 1, na.rm = TRUE) / 2
    maf_train_pred <- pmin(p_train_pred, 1 - p_train_pred)
    keep_markers <- maf_train_pred >= 0.01
    if (sum(keep_markers) < 100) keep_markers <- rep(TRUE, ncol(M_train))
    
    M_test_filtered <- M_test[, keep_markers, drop = FALSE]
    M_train_filtered <- M_train[, keep_markers, drop = FALSE]
    
    # Calculate cross-covariance K[test, train]
    # Using VanRaden formula components
    p <- colMeans(M_train_filtered + 1, na.rm = TRUE) / 2
    P <- 2 * (p - 0.5)
    denom <- sum(2 * p * (1 - p))
    if (denom < 1e-10) denom <- 1
    
    M_test_centered <- sweep(M_test_filtered, 2, P, "-")
    M_train_centered <- sweep(M_train_filtered, 2, P, "-")
    
    K_test_train <- (M_test_centered %*% t(M_train_centered)) / denom
    
    # Get K[train,train] from the fold
    K_train_train <- K_matrices_fold$A
    
    # Solve for prediction: u_test = K[test,train] * K[train,train]^-1 * u_train
    K_train_inv <- tryCatch({
      solve(K_train_train + diag(ridge, nrow(K_train_train)))
    }, error = function(e) {
      # Use pseudoinverse if singular
      MASS_ginv <- function(X) {
        svd_X <- svd(X)
        d_inv <- ifelse(svd_X$d > 1e-10, 1/svd_X$d, 0)
        svd_X$v %*% diag(d_inv) %*% t(svd_X$u)
      }
      MASS_ginv(K_train_train + diag(ridge, nrow(K_train_train)))
    })
    
    # Get training GEBVs as vector (average across environments)
    gebv_by_gid <- aggregate(GEBV ~ GID, data = gebv_train, FUN = mean, na.rm = TRUE)
    gebv_train_vec <- gebv_by_gid$GEBV[match(train_gids, gebv_by_gid$GID)]
    gebv_train_vec[is.na(gebv_train_vec)] <- 0
    
    # Predict test GEBVs
    gebv_test_vec <- as.vector(K_test_train %*% K_train_inv %*% gebv_train_vec)
    names(gebv_test_vec) <- test_gids
    
    # Create GEBV data frame for test individuals
    gebv_test <- data.frame(
      GID = test_gids,
      GEBV = gebv_test_vec,
      stringsAsFactors = FALSE
    )
  } else {
    # Use model's built-in predictions for all individuals
    gebv_by_gid <- aggregate(GEBV ~ GID, data = gebv_train, FUN = mean, na.rm = TRUE)
    gebv_test <- gebv_by_gid[gebv_by_gid$GID %in% test_gids, , drop = FALSE]
  }
  
  # Get test set phenotypes
  pheno_test <- pheno_full[fold$test, , drop = FALSE]
  pheno_test$GID <- as.character(pheno_test$GID)
  pheno_test$ENV <- as.character(pheno_test$ENV)
  
  # Merge predictions with test phenotypes
  test_with_pred <- merge(
    pheno_test[, c("GID", "ENV", "trait")],
    gebv_test,
    by = "GID",
    all.x = TRUE
  )
  
  # Add GxE breeding values if available
  # GxE effects (AE, DE, AAE) have BLUPs indexed by GID.ENV
  # These should be added to the main genetic GEBV for proper prediction
  if (any(c("AE", "DE", "AAE") %in% effects)) {
    fit <- model@model
    
    # Determine effect indices based on the order in build_design_matrices
    # The order is: A, D, ENV, AE, DE, AA, AAE
    effect_idx <- 0
    
    for (eff in c("A", "D", "ENV", "AE", "DE", "AA", "AAE")) {
      if (eff %in% effects) {
        effect_idx <- effect_idx + 1
        
        # Extract GxE BLUPs for this effect
        if (eff %in% c("AE", "DE", "AAE")) {
          u_gxe <- fit$u[[effect_idx]]
          if (!is.null(u_gxe)) {
            # GxE BLUPs are indexed by GID.ENV
            gxe_names <- names(u_gxe)
            if (is.null(gxe_names)) gxe_names <- rownames(u_gxe)
            
            if (!is.null(gxe_names)) {
              # Create GID.ENV identifiers for test observations
              test_gid_env <- paste(test_with_pred$GID, test_with_pred$ENV, sep = ".")
              
              # Match test observations to GxE levels
              gxe_match <- match(test_gid_env, gxe_names)
              valid_match <- !is.na(gxe_match)
              
              if (any(valid_match)) {
                gxe_vals <- as.numeric(u_gxe[gxe_match[valid_match]])
                test_with_pred$GEBV[valid_match] <- test_with_pred$GEBV[valid_match] + gxe_vals
              }
            }
          }
        }
      }
    }
  }
  
  # Extract matched predictions
  observed <- test_with_pred$trait
  predicted <- test_with_pred$GEBV
  test_ids <- test_with_pred$GID
  
  # Check for valid predictions
  valid <- !is.na(observed) & !is.na(predicted)
  
  if (sum(valid) == 0) {
    stop("No valid predictions for test set. Check GID matching.\n",
         "  Test GIDs sample: ", paste(head(unique(pheno_test$GID), 5), collapse = ", "), "\n",
         "  Predicted GIDs sample: ", paste(head(unique(gebv_test$GID), 5), collapse = ", "))
  }
  
  observed <- observed[valid]
  predicted <- predicted[valid]
  test_ids <- test_ids[valid]
  
  # Calculate predictive ability (correlation) - THE PRIMARY METRIC
  # This is cor(GEBV, y) - the standard metric in genomic prediction literature
  predictive_ability <- calc_accuracy(observed, predicted)
  mse <- calc_mse(observed, predicted)
  
  # Note: We do NOT calculate "true accuracy" here because:
  # 1. h² estimated from training data has high variance in small samples
  # 2. The formula accuracy = r/sqrt(h²) assumes TRUE h², not an estimate
  # 3. This would produce misleading/inflated values
  # 
  # Users who know the true h² (e.g., from simulation) can use:
  #   calc_prediction_accuracy(predictive_ability, h2_true)
  
  return(list(
    observed = observed,
    predicted = predicted,
    test_ids = test_ids,
    predictive_ability = predictive_ability,  # cor(GEBV, y) - PRIMARY METRIC
    mse = mse,
    n_test = length(observed),
    n_train = length(train_gids),
    recompute_K = recompute_K,
    model = model
  ))
}

#' Summarize Cross-Validation Results
#'
#' @param all_results List of results from all replications
#' @param scheme CV scheme used
#'
#' @return List with summary statistics
#' @keywords internal
summarize_cv_results <- function(all_results, scheme) {
  
  # Flatten results
  flat_results <- unlist(all_results, recursive = FALSE)
  
  # Remove failed folds
  valid_results <- flat_results[sapply(flat_results, function(x) {
    !is.null(x$predictive_ability) && !is.na(x$predictive_ability)
  })]
  
  if (length(valid_results) == 0) {
    return(list(
      mean_predictive_ability = NA,
      sd_predictive_ability = NA,
      mean_mse = NA,
      sd_mse = NA,
      metrics = data.frame()
    ))
  }
  
  # Extract metrics - predictive ability is THE primary metric
  predictive_abilities <- sapply(valid_results, function(x) x$predictive_ability)
  mses <- sapply(valid_results, function(x) x$mse)
  
  # Calculate summary statistics
  metrics <- data.frame(
    fold = sapply(valid_results, function(x) x$fold),
    replication = sapply(valid_results, function(x) x$replication),
    predictive_ability = predictive_abilities,
    mse = mses,
    n_test = sapply(valid_results, function(x) x$n_test),
    stringsAsFactors = FALSE
  )
  
  return(list(
    mean_predictive_ability = mean(predictive_abilities, na.rm = TRUE),
    sd_predictive_ability = sd(predictive_abilities, na.rm = TRUE),
    mean_mse = mean(mses, na.rm = TRUE),
    sd_mse = sd(mses, na.rm = TRUE),
    metrics = metrics
  ))
}

#' Compile Predictions from All Folds
#'
#' @param all_results List of results from all replications
#'
#' @return Data frame with all predictions
#' @keywords internal
compile_predictions <- function(all_results) {
  
  # Flatten results
  flat_results <- unlist(all_results, recursive = FALSE)
  
  # Remove failed folds
  valid_results <- flat_results[sapply(flat_results, function(x) {
    !is.null(x$observed) && length(x$observed) > 0
  })]
  
  if (length(valid_results) == 0) {
    return(data.frame(
      GID = character(),
      observed = numeric(),
      predicted = numeric(),
      fold = integer(),
      replication = integer(),
      stringsAsFactors = FALSE
    ))
  }
  
  # Compile into data frame
  predictions <- do.call(rbind, lapply(valid_results, function(x) {
    data.frame(
      GID = x$test_ids,
      observed = x$observed,
      predicted = x$predicted,
      fold = x$fold,
      replication = x$replication,
      stringsAsFactors = FALSE
    )
  }))
  
  return(predictions)
}

#' Plot Cross-Validation Results
#'
#' @param cv_result CVResult object
#' @param type Type of plot: "scatter", "boxplot", or "density"
#' @param ... Additional arguments passed to plotting functions
#'
#' @return Invisible NULL
#' @export
#' @examples
#' \dontrun{
#' cv_res <- cv_gblup(gdata, K, effects = "A", scheme = "CV1")
#' plot_cv_results(cv_res, type = "scatter")
#' plot_cv_results(cv_res, type = "boxplot")
#' }
plot_cv_results <- function(cv_result, type = c("scatter", "boxplot", "density"), ...) {
  
  type <- match.arg(type)
  
  if (!inherits(cv_result, "CVResult")) {
    stop("cv_result must be a CVResult object")
  }
  
  predictions <- cv_result@predictions
  
  if (nrow(predictions) == 0) {
    stop("No predictions available to plot")
  }
  
  if (type == "scatter") {
    # Scatter plot of observed vs predicted
    plot(predictions$observed, predictions$predicted,
         xlab = "Observed", ylab = "Predicted",
         main = paste("Cross-Validation Results -", cv_result@scheme),
         pch = 16, col = rgb(0, 0, 1, 0.3),
         ...)
    abline(0, 1, col = "red", lwd = 2)
    
    # Add correlation
    cor_val <- cor(predictions$observed, predictions$predicted, 
                   use = "complete.obs")
    legend("topleft", 
           legend = paste("r =", round(cor_val, 3)),
           bty = "n")
    
  } else if (type == "boxplot") {
    # Boxplot of predictive ability by fold
    metrics <- cv_result@metrics
    boxplot(predictive_ability ~ fold, data = metrics,
            xlab = "Fold", ylab = "Predictive Ability",
            main = paste("Predictive Ability by Fold -", cv_result@scheme),
            ...)
    
  } else if (type == "density") {
    # Density plot of residuals
    residuals <- predictions$observed - predictions$predicted
    plot(density(residuals, na.rm = TRUE),
         main = paste("Residual Distribution -", cv_result@scheme),
         xlab = "Residuals",
         ...)
    abline(v = 0, col = "red", lty = 2)
  }
  
  invisible(NULL)
}

#' Extract Cross-Validation Metrics
#'
#' @param cv_result CVResult object
#' @param by_fold Logical, return metrics by fold
#'
#' @return Data frame of metrics
#' @export
#' @examples
#' \dontrun{
#' cv_res <- cv_gblup(gdata, K, effects = "A", scheme = "CV1")
#' get_cv_metrics(cv_res)
#' get_cv_metrics(cv_res, by_fold = TRUE)
#' }
get_cv_metrics <- function(cv_result, by_fold = FALSE) {
  
  if (!inherits(cv_result, "CVResult")) {
    stop("cv_result must be a CVResult object")
  }
  
  if (by_fold) {
    return(cv_result@metrics)
  } else {
    # Overall summary
    metrics <- cv_result@metrics
    
    summary_df <- data.frame(
      scheme = cv_result@scheme,
      n_folds = cv_result@n_folds,
      n_reps = cv_result@n_reps,
      mean_predictive_ability = mean(metrics$predictive_ability, na.rm = TRUE),
      sd_predictive_ability = sd(metrics$predictive_ability, na.rm = TRUE),
      mean_mse = mean(metrics$mse, na.rm = TRUE),
      sd_mse = sd(metrics$mse, na.rm = TRUE),
      total_observations = sum(metrics$n_test),
      stringsAsFactors = FALSE
    )
    
    return(summary_df)
  }
}

#' Compare Multiple Models via Cross-Validation
#'
#' Runs cross-validation for multiple effect combinations and returns a comparison
#' table with statistical tests for significant differences.
#'
#' @param gdata GBLUPData object
#' @param K_matrices List of relationship matrices
#' @param effect_sets List of effect vectors to compare
#' @param scheme CV scheme ("CV1", "CV2", "CV0")
#' @param n_folds Number of folds (CV0 and CV1)
#' @param n_reps Number of replications
#' @param verbose Print progress
#'
#' @return Data frame comparing models
#' @importFrom stats t.test
#' @export
#' @examples
#' \dontrun{
#' K <- calc_relationship_matrices(gdata, matrices = c("A", "D"))
#' effects_to_compare <- list(
#'   "Additive" = "A",
#'   "Add+Dom" = c("A", "D"),
#'   "Add+Dom+Epi" = c("A", "D", "AA")
#' )
#' comparison <- compare_models_cv(gdata, K, effects_to_compare)
#' print(comparison)
#' }
compare_models_cv <- function(gdata, K_matrices, effect_sets,
                               scheme = "CV1", n_folds = 5, n_reps = 1,
                               verbose = TRUE) {
  
  if (!is.list(effect_sets)) {
    stop("effect_sets must be a named list of effect vectors")
  }
  
  model_names <- names(effect_sets)
  if (is.null(model_names)) {
    model_names <- paste0("Model_", seq_along(effect_sets))
  }
  
  results <- list()
  
  for (i in seq_along(effect_sets)) {
    model_name <- model_names[i]
    effects <- effect_sets[[i]]
    
    if (verbose) {
      cat("\n=== Evaluating:", model_name, "===\n")
      cat("Effects:", paste(effects, collapse = ", "), "\n")
    }
    
    cv_result <- tryCatch({
      cv_gblup(gdata, K_matrices, effects = effects, 
               scheme = scheme, n_folds = n_folds, n_reps = n_reps,
               verbose = verbose)
    }, error = function(e) {
      warning("Model '", model_name, "' failed: ", e$message)
      NULL
    })
    
    if (!is.null(cv_result)) {
      metrics <- cv_result@metrics
      results[[model_name]] <- list(
        effects = paste(effects, collapse = "+"),
        mean_r = mean(metrics$predictive_ability, na.rm = TRUE),
        sd_r = sd(metrics$predictive_ability, na.rm = TRUE),
        mean_mse = mean(metrics$mse, na.rm = TRUE),
        sd_mse = sd(metrics$mse, na.rm = TRUE),
        n_folds = nrow(metrics),
        fold_correlations = metrics$predictive_ability
      )
    }
  }
  
  if (length(results) == 0) {
    stop("All models failed")
  }
  
  # Create comparison table
  comparison <- data.frame(
    model = names(results),
    effects = sapply(results, `[[`, "effects"),
    mean_r = sapply(results, `[[`, "mean_r"),
    sd_r = sapply(results, `[[`, "sd_r"),
    mean_mse = sapply(results, `[[`, "mean_mse"),
    sd_mse = sapply(results, `[[`, "sd_mse"),
    n_folds = sapply(results, `[[`, "n_folds"),
    stringsAsFactors = FALSE
  )
  
  # Rank by mean predictive ability
  comparison <- comparison[order(-comparison$mean_r), ]
  comparison$rank <- seq_len(nrow(comparison))
  
  # Add pairwise t-tests against best model if more than one model
  if (nrow(comparison) > 1) {
    best_model <- comparison$model[1]
    best_r <- results[[best_model]]$fold_correlations
    
    comparison$p_vs_best <- sapply(comparison$model, function(m) {
      if (m == best_model) return(NA)
      m_r <- results[[m]]$fold_correlations
      if (length(m_r) < 2 || length(best_r) < 2) return(NA)
      tryCatch({
        t.test(best_r, m_r, paired = TRUE)$p.value
      }, error = function(e) NA)
    })
    
    comparison$sig <- ifelse(is.na(comparison$p_vs_best), "",
                              ifelse(comparison$p_vs_best < 0.001, "***",
                                     ifelse(comparison$p_vs_best < 0.01, "**",
                                            ifelse(comparison$p_vs_best < 0.05, "*", ""))))
  }
  
  if (verbose) {
    cat("\n=== Model Comparison ===\n")
    print(comparison[, c("rank", "model", "effects", "mean_r", "sd_r", "sig")])
    if (nrow(comparison) > 1) {
      cat("\nSignificance codes: *** p<0.001, ** p<0.01, * p<0.05\n")
      cat("p-values from paired t-test vs. best model\n")
    }
  }
  
  return(comparison)
}