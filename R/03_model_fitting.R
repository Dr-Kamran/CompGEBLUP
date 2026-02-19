#' Fit GBLUP Model
#'
#' Fits a Genomic Best Linear Unbiased Prediction (GBLUP) model
#' using the internal REML-based mixed model engine with optional
#' EMMA-style efficient algorithm for single random effect models.
#'
#' @param gdata GBLUPData object
#' @param K_matrices List of relationship matrices
#' @param effects Character vector of random effects to include
#' @param ridge Ridge parameter for numerical stability (default: 0.001).
#'   Applied uniformly to all K matrices unless \code{ridge_per_matrix} is provided.
#' @param ridge_per_matrix Optional named numeric vector of per-matrix ridge values.
#'   When provided, each K matrix receives its designated ridge value instead of the
#'   scalar \code{ridge}. Names must match effect names (e.g. "A", "D", "AA", "AE").
#'   Effects not listed fall back to the scalar \code{ridge}.
#'   Example: \code{ridge_per_matrix = c(A = 0.001, D = 0.05, AA = 0.05)}.
#' @param nIters Maximum iterations for REML (default: 100). Increase for complex
#'   models with 3+ random effects if convergence warnings occur.
#' @param tolParConvLL Convergence tolerance
#' @param verbose Logical, print progress messages
#' @param use.emma Control solver selection (default TRUE).
#'   When TRUE: auto-selects the fastest solver (EMMA for single-effect small n,
#'   Henderson EM-REML for large n >> q). When FALSE: uses Newton-based solvers
#'   (Henderson AI-REML for large n >> q, observation-space AI-REML for small n).
#'   Set to FALSE if variance estimates drift or fail to converge with the default.
#' @param use.sparse Use sparse matrices for design matrices (default FALSE)
#' @param K_already_regularized Logical, if TRUE skips adding ridge to K matrices.
#'   Set to TRUE if K matrices were already regularized via make_positive_definite()
#'   to avoid double regularization.
#' @param fixed Optional fixed effects specification. Can be NULL (intercept only, default),
#'   a one-sided formula (e.g., ~ TESTER or ~ TESTER + BLOCK) referencing columns in 
#'   gdata@@phenotypes, or a pre-built design matrix. When provided, fixed effects are
#'   estimated jointly with random effects in the mixed model.
#'
#' @return GBLUPModel object
#' @export
#' @examples
#' \dontrun{
#' gdata <- simulate_genotypes(n_ind = 100, n_snp = 500)
#' gdata <- simulate_phenotypes(gdata, n_env = 3, h2 = 0.6)
#' K <- calc_relationship_matrices(gdata, matrices = "A")
#' model <- fit_gblup(gdata, K, effects = "A")
#' summary(model)
#' 
#' # With differential ridge
#' model_diff <- fit_gblup(gdata, K, effects = c("A", "D"),
#'                         ridge_per_matrix = c(A = 0.001, D = 0.05))
#' }
fit_gblup <- function(gdata, K_matrices, effects = "A", 
                      ridge = 0.001,
                      ridge_per_matrix = NULL,
                      nIters = 100,
                      tolParConvLL = 1e-4, verbose = TRUE,
                      use.emma = TRUE, use.sparse = FALSE,
                      K_already_regularized = FALSE,
                      fixed = NULL) {
  
  # ===== INPUT VALIDATION =====
  if (!inherits(gdata, "GBLUPData")) {
    stop("gdata must be a GBLUPData object")
  }
  
  if (!is.list(K_matrices) || length(K_matrices) == 0) {
    stop("K_matrices must be a non-empty list of relationship matrices")
  }
  
  if (!is.character(effects) || length(effects) == 0) {
    stop("effects must be a non-empty character vector")
  }
  
  # Validate effects
  builtin_effects <- c("A", "D", "ENV", "AE", "DE", "AA", "AAE")
  custom_effects <- setdiff(effects, builtin_effects)
  if (length(custom_effects) > 0) {
    # Custom effects are allowed if they have a matching kernel in K_matrices
    missing_kernels <- setdiff(custom_effects, names(K_matrices))
    if (length(missing_kernels) > 0) {
      stop("Custom effects without matching K_matrices: ",
           paste(missing_kernels, collapse = ", "),
           ". Either provide kernels in K_matrices or use built-in effects: ",
           paste(builtin_effects, collapse = ", "))
    }
    if (verbose) {
      cat("Custom effects detected:", paste(custom_effects, collapse = ", "), "\n")
    }
  }
  
  # ===== PREPARE DATA =====
  pheno <- gdata@phenotypes
  
  # Ensure character types for factors
  pheno$GID <- as.character(pheno$GID)
  pheno$ENV <- as.character(pheno$ENV)
  
  # Store ALL unique GIDs from genotypes (for prediction)
  all_gids <- rownames(gdata@genotypes)
  
  # Keep track of which observations have phenotypes
  has_pheno <- !is.na(pheno$trait)
  
  # Get phenotyped subset for model fitting
  pheno_complete <- pheno[has_pheno, ]
  
  if (nrow(pheno_complete) < 10) {
    stop("Too few complete observations (", nrow(pheno_complete), ")")
  }
  
  # Get response (only complete cases)
  y <- as.matrix(pheno_complete$trait)
  
  # Check for signs of pre-corrected data
  trait_mean <- mean(pheno_complete$trait, na.rm = TRUE)
  trait_sd <- sd(pheno_complete$trait, na.rm = TRUE)
  if (abs(trait_mean) < trait_sd * 0.01 && trait_sd > 0) {
    message("Note: Trait mean is near zero (", round(trait_mean, 4), 
            "), which may indicate pre-corrected/residualized data. ",
            "If you used residuals from a fixed-effects model (e.g., lm(yield ~ Tester)), ",
            "consider using the 'fixed' parameter for joint estimation instead, ",
            "which is statistically superior. See ?fit_gblup for details.")
  }
  
  # Count environments for heritability calculation
  n_env <- length(unique(pheno_complete$ENV))
  
  # ===== BUILD DESIGN MATRICES =====
  # Pass un-ridged K matrices so Kronecker products are computed from clean matrices
  design <- build_design_matrices(pheno_complete, K_matrices, effects, 
                                   use.sparse = use.sparse, fixed = fixed)
  
  # Add ridge to K matrices AFTER design matrix construction
  # This ensures auto-computed GxE matrices also get proper ridge
  if (!K_already_regularized) {
    # BUG FIX (v2.3.1): Support ridge_per_matrix for differential regularisation.
    # Previously applied scalar ridge uniformly to all K matrices.
    .get_ridge_fit <- function(mat_type) {
      if (!is.null(ridge_per_matrix) && mat_type %in% names(ridge_per_matrix)) {
        ridge_per_matrix[[mat_type]]
      } else {
        ridge
      }
    }
    for (k_idx in seq_along(design$K_list)) {
      k_name <- design$effect_names[k_idx]
      r_val <- .get_ridge_fit(k_name)
      design$K_list[[k_idx]] <- design$K_list[[k_idx]] +
        diag(r_val, nrow(design$K_list[[k_idx]]), ncol(design$K_list[[k_idx]]))
    }
    if (verbose) {
      if (!is.null(ridge_per_matrix)) {
        cat("Differential ridge applied to K matrices:\n")
        for (k_idx in seq_along(design$K_list)) {
          k_name <- design$effect_names[k_idx]
          cat(sprintf("  %s: %.4f\n", k_name, .get_ridge_fit(k_name)))
        }
      } else {
        cat("Ridge parameter added to all K matrices:", ridge, "\n")
      }
    }
  } else {
    if (verbose) {
      cat("Using pre-regularized K matrices\n")
    }
  }
  
  X <- design$X
  Z_list <- design$Z_list
  K_list <- design$K_list
  
  # Check X matrix rank and remove collinear columns
  x_rank <- qr(X)$rank
  if (x_rank < ncol(X)) {
    pivot <- qr(X)$pivot[seq_len(x_rank)]
    dropped <- colnames(X)[setdiff(seq_len(ncol(X)), pivot)]
    X <- X[, pivot, drop = FALSE]
    warning("Fixed effects matrix is rank-deficient. Dropped collinear column(s): ",
            paste(dropped, collapse = ", "), call. = FALSE)
  }
  
  if (length(Z_list) == 0) {
    stop("No random effects specified. At least one effect is required.")
  }
  
  if (verbose) {
    cat("Fitting GBLUP model with effects:", paste(effects, collapse = ", "), "\n")
    cat("Number of observations:", nrow(y), "\n")
    cat("Number of random effects:", length(Z_list), "\n")
    for (i in seq_along(Z_list)) {
      cat("  ", design$effect_names[i], ":", ncol(Z_list[[i]]), "levels\n")
    }
  }
  
  # ===== FIT MODEL =====
  # BUG FIX (v2.3.1): tolParInv is for internal matrix inversion stability,
  # not for regularising K matrices. Use a small constant (1e-6) rather than
  # the user's ridge parameter which can be large (e.g. 0.05 for AE models).
  fit <- tryCatch({
    fit_mme(
      y = y,
      X = X,
      Z_list = Z_list,
      K_list = K_list,
      nIters = nIters,
      tolParConvLL = tolParConvLL,
      tolParInv = 1e-6,
      verbose = verbose,
      use.emma = use.emma
    )
  }, error = function(e) {
    stop("Model fitting failed: ", e$message, "\n",
         "  Check that:\n",
         "  1. Data has sufficient variation\n",
         "  2. Relationship matrices are valid\n",
         "  3. No singular matrices (try increasing ridge parameter)")
  })
  
  # ===== CHECK CONVERGENCE =====
  if (!fit$convergence && verbose) {
    warning("Model did not converge. Results may be unreliable.\n",
            "  Try increasing nIters or ridge parameter.",
            call. = FALSE)
  }
  
  # ===== EXTRACT VARIANCE COMPONENTS =====
  theta <- fit$theta
  effect_names_with_res <- c(design$effect_names, "Residual")
  
  varcomp <- data.frame(
    Component = effect_names_with_res,
    Variance = theta,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  # Check for negative variances
  if (any(theta < 0)) {
    warning("Negative variance components detected. ",
            "Consider increasing ridge parameter.",
            call. = FALSE)
  }
  
  # Check if any genetic variance component hit the solver boundary (near zero).
  # Two-tier detection:
  #
  # 1. Absolute threshold (all models): θ ≤ MIN_VARIANCE × 10 (= 1e-9).
  #    Fires when Newton-based solvers clamp to MIN_VARIANCE floor.
  #
  # 2. Proportional threshold (single genetic component only): θ/total < 0.5%.
  #    EM-REML has a positive bias for σ²_A near zero (the trace correction
  #    tr(C⁻¹K⁻¹)/q creates a floor at approximately σ²_A ≈ σ²_e × q/n).
  #    For single-component models (effects="A"), this proportional check
  #    detects when σ²_A converges to a negligible fraction of total variance.
  #    NOT applied to multi-component models because interaction terms (AE, DE)
  #    are legitimately small — warning on those would be a false positive.
  #
  # Threshold rationale: 0.5% is well below any biologically meaningful h².
  # All expect_silent tests use h2 ≥ 0.5 (100× above threshold).
  genetic_theta <- theta[seq_len(length(theta) - 1)]
  n_genetic <- length(genetic_theta)
  total_var <- sum(theta)
  
  at_boundary <- any(genetic_theta <= .MME_CONSTANTS$MIN_VARIANCE * 10)
  
  # Proportional check: only for single-component models
  if (!at_boundary && n_genetic == 1 && total_var > 0) {
    if (genetic_theta[1] / total_var < 5e-3) {
      at_boundary <- TRUE
    }
  }
  
  # Also check solver's own boundary tracking (AI-REML)
  if (!at_boundary && !is.null(fit$boundary)) {
    at_boundary <- any(fit$boundary[seq_len(n_genetic)])
  }
  
  if (at_boundary) {
    warning("One or more genetic variance components are at or near zero (boundary). ",
            "This often indicates:\n",
            "  1. Pre-corrected/residualized phenotypes that lost genetic signal\n",
            "  2. Model misspecification\n",
            "  3. No genetic variation in the trait\n",
            "If you pre-corrected phenotypes (e.g., using residuals from lm()), ",
            "consider using the 'fixed' parameter instead for joint estimation:\n",
            "  fit_gblup(..., fixed = ~ TESTER)\n",
            "See ?fit_gblup for details.",
            call. = FALSE)
  }
  
  if (verbose) {
    cat("\nVariance components:\n")
    print(varcomp)
    
    # Print heritability if A effect is present
    # Note: For multi-environment models, use heritability() function for proper estimates
    # that account for G×E interaction variance correctly
    if ("A" %in% design$effect_names) {
      var_a <- theta[which(design$effect_names == "A")]
      var_total <- sum(theta)
      h2_simple <- var_a / var_total
      cat("\nHeritability (h2, simplified estimate):", round(h2_simple, 3), "\n")
      if (n_env > 1) {
        cat("  Note: For multi-environment models, use heritability() function for\n")
        cat("  proper estimates that correctly handle G×E variance components.\n")
      }
    }
  }
  
  # ===== EXTRACT GEBVs FOR ALL INDIVIDUALS =====
  gebv_df <- extract_gebvs_from_fit(fit, pheno, design, effects, K_matrices, all_gids)
  
  if (verbose) {
    cat("\nExtracted", nrow(gebv_df), "GEBVs\n")
    cat("GEBV range:", round(range(gebv_df$GEBV, na.rm = TRUE), 3), "\n")
  }
  
  # Store n_env in fit for heritability calculation
  fit$n_env <- n_env
  
  # ===== CREATE MODEL OBJECT =====
  model <- new("GBLUPModel",
               model = fit,
               data = gdata,
               gebv = gebv_df,
               varcomp = varcomp,
               description = paste("GBLUP:", paste(effects, collapse = "+")),
               effects = effects
  )
  
  return(model)
}

#' Extract GEBVs from fitted model using proper BLUP equations
#'
#' @param fit Fitted model from fit_mme
#' @param pheno Phenotype data (full, including NAs)
#' @param design Design matrices from build_design_matrices
#' @param effects Effects included
#' @param K_matrices Un-regularized K matrices (for BLUP prediction)
#' @param all_gids All GIDs from genotype matrix
#' @return Data frame with GEBVs
#' @keywords internal
extract_gebvs_from_fit <- function(fit, pheno, design, effects, K_matrices = NULL, all_gids = NULL) {
  
  # Find the main genetic effect (A, D, AA, or custom GID-level kernels)
  builtin_genetic <- c("A", "D", "AA")
  non_genetic <- c("ENV", "AE", "DE", "AAE")
  
  genetic_idx <- which(design$effect_names %in% builtin_genetic)[1]
  
  if (is.na(genetic_idx)) {
    # Look for custom GID-level effects (not ENV or interaction effects)
    custom_candidates <- which(!(design$effect_names %in% non_genetic))
    if (length(custom_candidates) > 0) {
      genetic_idx <- custom_candidates[1]
    } else {
      warning("Could not identify genetic effect. Using first random effect.")
      genetic_idx <- 1
    }
  }
  
  # Get BLUPs for genetic effect
  u_genetic <- fit$u[[genetic_idx]]
  
  # Get the corresponding Z matrix columns (GID names) - these are the TRAINED GIDs
  Z_genetic <- design$Z_list[[genetic_idx]]
  trained_gids <- as.character(colnames(Z_genetic))
  
  # Ensure u_genetic is a vector
  if (is.matrix(u_genetic)) {
    u_genetic <- as.vector(u_genetic)
  }
  
  # Create GEBV lookup for trained individuals
  names(u_genetic) <- trained_gids
  
  gebv_trained <- data.frame(
    GID = trained_gids,
    GEBV = as.numeric(u_genetic),
    stringsAsFactors = FALSE
  )
  
  # Predict GEBVs for ALL individuals using proper BLUP equations
  # Find the kernel matching the genetic effect
  genetic_effect_name <- design$effect_names[genetic_idx]
  genetic_kernel_name <- if (genetic_effect_name %in% names(K_matrices)) {
    genetic_effect_name
  } else if ("A" %in% names(K_matrices)) {
    "A"
  } else {
    NULL
  }
  
  if (!is.null(K_matrices) && !is.null(all_gids) && !is.null(genetic_kernel_name)) {
    K_A <- K_matrices[[genetic_kernel_name]]
    all_gids <- as.character(all_gids)
    
    # Find individuals NOT in training set but in K matrix
    untrained_gids <- setdiff(all_gids, trained_gids)
    
    if (length(untrained_gids) > 0) {
      # Check that GIDs exist in K matrix
      valid_untrained <- untrained_gids[untrained_gids %in% rownames(K_A)]
      valid_trained <- trained_gids[trained_gids %in% rownames(K_A)]
      
      if (length(valid_untrained) > 0 && length(valid_trained) > 0) {
        # Use proper BLUP prediction function
        gebv_train_vec <- u_genetic[valid_trained]
        
        gebv_untrained_vals <- predict_blup_new(
          u_train = gebv_train_vec,
          K = K_A,
          train_ids = valid_trained,
          new_ids = valid_untrained
        )
        
        gebv_untrained <- data.frame(
          GID = valid_untrained,
          GEBV = gebv_untrained_vals,
          stringsAsFactors = FALSE
        )
        
        # Combine trained and untrained
        gebv_all <- rbind(gebv_trained, gebv_untrained)
      } else {
        gebv_all <- gebv_trained
      }
    } else {
      gebv_all <- gebv_trained
    }
  } else {
    gebv_all <- gebv_trained
  }
  
  # Ensure GID is character
  gebv_all$GID <- as.character(gebv_all$GID)
  
  # Prepare phenotype data
  pheno_merge <- data.frame(
    GID = as.character(pheno$GID),
    ENV = as.character(pheno$ENV),
    observed = pheno$trait,
    stringsAsFactors = FALSE
  )
  
  # Merge with phenotype data
  gebv_df <- merge(pheno_merge, gebv_all, by = "GID", all.x = TRUE)
  
  # Reorder columns
  gebv_df <- gebv_df[, c("GID", "ENV", "GEBV", "observed")]
  
  return(gebv_df)
}

#' Predict breeding values for new individuals
#'
#' Uses proper BLUP equations for prediction of new individuals.
#'
#' @param model GBLUPModel object
#' @param newdata Optional new phenotype data
#' @param K_new Optional relationship matrix for new individuals (rows = new, cols = training)
#'
#' @return Data frame with predicted breeding values
#' @export
#' @examples
#' \dontrun{
#' model <- fit_gblup(gdata, K, effects = "A")
#' predictions <- predict_gblup(model)
#' }
predict_gblup <- function(model, newdata = NULL, K_new = NULL) {
  
  if (!inherits(model, "GBLUPModel")) {
    stop("model must be a GBLUPModel object")
  }
  
  if (is.null(newdata) && is.null(K_new)) {
    return(model@gebv[, c("GID", "ENV", "GEBV")])
  }
  
  # For new data prediction with K_new
  if (!is.null(K_new)) {
    # Get training GEBVs
    gebv_train <- model@gebv
    gebv_by_gid <- aggregate(GEBV ~ GID, data = gebv_train, FUN = mean, na.rm = TRUE)
    
    gid_train <- as.character(gebv_by_gid$GID)
    gid_new <- rownames(K_new)
    if (is.null(gid_new)) {
      gid_new <- paste0("NEW", seq_len(nrow(K_new)))
    }
    
    # Match columns with training individuals
    common_train <- intersect(colnames(K_new), gid_train)
    if (length(common_train) == 0) {
      stop("No overlap between new K matrix columns and training GIDs")
    }
    
    # Get original K matrix for training individuals
    K_orig <- calc_relationship_matrices(model@data, matrices = "A")$A
    
    # Use proper BLUP prediction: u_new = K_new,train * K_train,train^{-1} * u_train
    # We only need K_new,train (the cross-relationship) and K_train,train
    gebv_train_common <- gebv_by_gid$GEBV[match(common_train, gebv_by_gid$GID)]
    
    # Construct a minimal K matrix with only the blocks needed for prediction
    # Note: predict_blup_new will extract K[new_ids, train_ids] and K[train_ids, train_ids]
    # We don't need K[new_ids, new_ids] for prediction
    combined_size <- length(gid_new) + length(common_train)
    combined_K <- matrix(0, nrow = combined_size, ncol = combined_size)
    rownames(combined_K) <- c(gid_new, common_train)
    colnames(combined_K) <- c(gid_new, common_train)
    
    # Fill in the blocks we have:
    # K[train, train] - from original K
    combined_K[common_train, common_train] <- K_orig[common_train, common_train]
    # K[new, train] - from K_new (cross-relationship)
    combined_K[gid_new, common_train] <- K_new[, common_train]
    # K[train, new] - symmetric
    combined_K[common_train, gid_new] <- t(K_new[, common_train])
    # K[new, new] - left as zeros (not used by predict_blup_new)
    
    gebv_new <- predict_blup_new(
      u_train = gebv_train_common,
      K = combined_K,
      train_ids = common_train,
      new_ids = gid_new
    )
    
    result <- data.frame(
      GID = gid_new,
      ENV = NA_character_,
      GEBV = gebv_new,
      stringsAsFactors = FALSE
    )
    
    return(result)
  }
  
  stop("Prediction for new individuals requires K_new matrix.")
}

#' Calculate heritability from GBLUP model
#'
#' Calculates narrow-sense (h²) and broad-sense (H²) heritability from 
#' variance components estimated by GBLUP.
#'
#' @details
#' For multi-environment models, heritability can be computed two ways:
#' 
#' **For selection on mean performance (default output):**
#' G×E interactions are treated as NOISE - they reduce selection accuracy.
#' \itemize{
#'   \item h2_mean = Va / (Va + Vae/n_env + Ve)
#'   \item H2_mean = (Va + Vd) / (Va + Vd + Vae/n_env + Vde/n_env + Ve)
#' }
#' Use this when selecting genotypes for stable performance across environments.
#' 
#' **For variance partitioning:**
#' G×E interactions counted as genetic variance.
#' \itemize{
#'   \item h2_total = (Va + Vae/n_env) / (Va + Vae/n_env + Ve)
#'   \item H2_total = (Va + Vd + Vae/n_env + Vde/n_env) / Vp
#' }
#' Use this for understanding genetic architecture and variance components.
#'
#' For single-environment models (n_env = 1), both give identical results.
#' 
#' @section Warning - Interpretation:
#' For multi-environment data, the default (h2, H2) assumes selection on
#' MEAN performance across environments. For single-environment selection
#' or variance partitioning, use h2_total and H2_total. When in doubt,
#' consult a quantitative geneticist.
#'
#' @param model GBLUPModel object
#' @param n_env Number of environments (auto-detected if NULL)
#' @return Named numeric vector with heritability estimates. The vector contains:
#'   \item{h2}{Narrow-sense heritability for mean performance (G×E as noise)}
#'   \item{H2}{Broad-sense heritability for mean performance (G×E as noise)}
#'   \item{h2_total}{Narrow-sense including G×E as genetic}
#'   \item{H2_total}{Broad-sense including G×E as genetic}
#'   
#'   The number of environments used in calculation is stored as an attribute 
#'   and can be accessed via \code{attr(result, "n_env")}.
#' @export
#' @examples
#' \dontrun{
#' model <- fit_gblup(gdata, K, effects = c("A", "AE"))
#' h <- heritability(model)
#' 
#' # For selection decisions:
#' h["h2"]       # Heritability of mean performance
#' 
#' # For variance partitioning:
#' h["h2_total"] # Total genetic variance explained
#' 
#' # Access number of environments:
#' attr(h, "n_env")
#' }
heritability <- function(model, n_env = NULL) {
  
  if (!inherits(model, "GBLUPModel")) {
    stop("model must be a GBLUPModel object")
  }
  
  vc <- model@varcomp
  
  # Find residual variance - should be exactly one
  res_idx <- grep("Residual", vc$Component, ignore.case = TRUE)
  if (length(res_idx) == 0) {
    warning("Could not identify residual variance component")
    result <- c(h2 = NA, H2 = NA)
    attr(result, "n_env") <- NA
    return(result)
  }
  if (length(res_idx) > 1) {
    warning("Multiple residual components found. Using sum.")
  }
  
  var_e <- sum(vc$Variance[res_idx])
  
  # Auto-detect n_env from model if not provided
  if (is.null(n_env)) {
    if (!is.null(model@model$n_env)) {
      n_env <- model@model$n_env
    } else if ("ENV" %in% names(model@gebv)) {
      n_env <- length(unique(model@gebv$ENV))
    } else {
      n_env <- 1
    }
  }
  
  # Identify variance component types (excluding residual)
  components <- vc$Component[-res_idx]
  variances <- vc$Variance[-res_idx]
  
  if (length(components) == 0) {
    warning("No genetic variance components found")
    result <- c(h2 = NA, H2 = NA)
    attr(result, "n_env") <- n_env
    return(result)
  }
  
  # Classify effects:
  # - Main genetic effects: A, D, AA (no E suffix)
  # - Interaction effects: AE, DE, AAE (end with E but not ENV or E)
  # - Environment effect: ENV or E (non-genetic variance)
  is_env <- components == "ENV" | components == "E"
  is_interaction <- grepl("E$", components) & !is_env
  is_main_genetic <- !is_interaction & !is_env
  
  # Main genetic variances
  var_main <- sum(variances[is_main_genetic])
  
  # Interaction variances (divided by n_env for contribution to mean)
  var_interaction <- sum(variances[is_interaction])
  var_interaction_contrib <- if (n_env > 1) var_interaction / n_env else var_interaction
  
  # Environment variance (non-genetic, contributes to phenotypic variance)
  var_env <- sum(variances[is_env])
  
  # Total phenotypic variance for mean across environments
  # NOTE: σ²_ENV is excluded from the denominator for entry-mean heritability.
  # When selecting on the genotype mean across environments, the environmental
  # main effect cancels out (each genotype is tested in all environments).
  # The relevant denominator is: σ²_A + σ²_GxE/n_env + σ²_e/(n_env * n_rep)
  # For unbalanced data, we approximate n_rep = 1 (one BLUE per GID x ENV).
  total_var_entry_mean <- var_main + var_interaction_contrib + var_e
  
  # Total phenotypic variance (all components, for variance partitioning)
  total_var_full <- var_main + var_interaction_contrib + var_env + var_e
  
  if (total_var_entry_mean <= 0) {
    warning("Total variance is zero or negative")
    result <- c(h2 = NA, H2 = NA, h2_total = NA, H2_total = NA)
    attr(result, "n_env") <- n_env
    return(result)
  }
  
  # === HERITABILITY FOR SELECTION ON MEAN (G×E as noise) ===
  # Entry-mean h²: σ²_ENV excluded from denominator
  H2_mean <- var_main / total_var_entry_mean
  
  # === HERITABILITY FOR VARIANCE PARTITIONING (G×E as genetic) ===
  # Uses full denominator including σ²_ENV for proportion-of-total
  H2_total <- (var_main + var_interaction_contrib) / total_var_full
  
  # Narrow-sense heritability
  a_idx <- which(components == "A")
  ae_idx <- which(components == "AE")
  
  if (length(a_idx) > 0) {
    var_a <- sum(variances[a_idx])
    var_ae <- if (length(ae_idx) > 0) sum(variances[ae_idx]) else 0
    var_ae_contrib <- if (n_env > 1 && var_ae > 0) var_ae / n_env else var_ae
    
    # For selection on mean: only main additive effect, entry-mean denominator
    h2_mean <- var_a / total_var_entry_mean
    
    # For variance partitioning: include A×E, full denominator
    h2_total <- (var_a + var_ae_contrib) / total_var_full
  } else {
    h2_mean <- NA
    h2_total <- NA
  }
  
  # Bound check (can exceed 1 due to estimation error)
  if (!is.na(h2_mean) && h2_mean > 1) {
    warning("h2 > 1 (", round(h2_mean, 3), "). Estimation error or model misspecification.")
  }
  if (!is.na(H2_mean) && H2_mean > 1) {
    warning("H2 > 1 (", round(H2_mean, 3), "). Estimation error or model misspecification.")
  }
  
  result <- c(
    h2 = h2_mean,           # Default: for selection on mean
    H2 = H2_mean,           # Default: for selection on mean
    h2_total = h2_total,    # For variance partitioning
    H2_total = H2_total     # For variance partitioning
  )
  attr(result, "n_env") <- n_env
  return(result)
}

#' Get variance components from model
#'
#' @param model GBLUPModel object
#' @return Data frame of variance components with columns Component and Variance
#' @export
#' @examples
#' \dontrun{
#' model <- fit_gblup(gdata, K, effects = "A")
#' vc <- get_varcomp(model)
#' print(vc)
#' }
get_varcomp <- function(model) {
  if (!inherits(model, "GBLUPModel")) {
    stop("model must be a GBLUPModel object")
  }
  return(model@varcomp)
}

#' Compute Prediction Error Variance (PEV)
#'
#' Computes an approximation of Prediction Error Variance for each individual.
#' 
#' @details
#' True PEV requires the full coefficient matrix inverse from the mixed model equations,
#' which is computationally expensive to store and compute for large datasets. This function
#' provides a practical approximation based on individual prediction quality:
#' 
#' PEV_i approximately equals sigma^2_g * (1 - r^2_i)
#' 
#' where r^2_i is the squared correlation between GEBV and observed phenotype for individual i
#' within their environment group. This approximation correctly identifies individuals with
#' higher/lower prediction reliability, though it may not match MME-based PEV exactly.
#' 
#' For small datasets (<1000 individuals), consider computing true PEV using the full
#' MME coefficient matrix inverse.
#'
#' @param model GBLUPModel object
#' @return Vector of PEV values for each individual
#' @export
compute_pev <- function(model) {
  
  if (!inherits(model, "GBLUPModel")) {
    stop("model must be a GBLUPModel object")
  }
  
  # Get genetic variance
  vc <- model@varcomp
  a_idx <- grep("^A$", vc$Component)
  
  if (length(a_idx) == 0) {
    warning("Additive variance component not found, using total genetic variance")
    # Use sum of all non-residual variance components
    res_idx <- which(vc$Component == "Residual")
    if (length(res_idx) > 0) {
      var_g <- sum(vc$Variance[-res_idx])
    } else {
      return(rep(NA, nrow(model@gebv)))
    }
  } else {
    var_g <- vc$Variance[a_idx]
  }
  
  if (var_g <= 0) {
    warning("Genetic variance is zero or negative")
    return(rep(NA, nrow(model@gebv)))
  }
  
  gebv <- model@gebv
  pev <- numeric(nrow(gebv))
  
  # Compute PEV for each environment group separately
  # This accounts for different amounts of information per environment
  env_groups <- unique(gebv$ENV)
  
  for (env in env_groups) {
    idx <- which(gebv$ENV == env)
    env_data <- gebv[idx, ]
    
    # Find individuals with both GEBV and observed values
    complete <- complete.cases(env_data$GEBV, env_data$observed)
    
    if (sum(complete) >= 2) {
      # Compute correlation between GEBV and observed for this environment
      r <- cor(env_data$GEBV[complete], env_data$observed[complete], use = "complete.obs")
      r_squared <- r^2
      
      # PEV approximation: higher correlation = lower PEV
      # Individuals with data: PEV = var_g * (1 - r²)
      pev[idx[complete]] <- var_g * (1 - r_squared)
      
      # Individuals without observed data in this environment get maximum PEV
      pev[idx[!complete]] <- var_g
    } else {
      # Not enough data for correlation - use maximum PEV
      pev[idx] <- var_g
    }
  }
  
  return(pev)
}

#' Compute Reliability of GEBVs
#'
#' Computes reliability (1 - PEV/σ²_g) for each individual's GEBV.
#' Reliability indicates the proportion of genetic variance captured by the prediction.
#' 
#' @details
#' Reliability ranges from 0 (no information) to 1 (perfect information). 
#' Values typically range from 0.3-0.8 depending on:
#' - Training population size
#' - Genetic relationship to training population  
#' - Heritability of the trait
#' - Marker density
#' 
#' This function uses an approximation of PEV based on prediction quality within
#' environment groups. For true MME-based reliability, the full coefficient matrix
#' inverse would be needed, which is computationally expensive for large datasets.
#'
#' @param model GBLUPModel object
#' @return Vector of reliability values (0 to 1)
#' @export
compute_reliability <- function(model) {
  
  if (!inherits(model, "GBLUPModel")) {
    stop("model must be a GBLUPModel object")
  }
  
  # Get genetic variance
  vc <- model@varcomp
  a_idx <- grep("^A$", vc$Component)
  
  if (length(a_idx) == 0) {
    warning("Additive variance component not found")
    return(rep(NA, nrow(model@gebv)))
  }
  
  var_a <- vc$Variance[a_idx]
  
  if (var_a <= 0) {
    warning("Additive variance is zero or negative")
    return(rep(NA, nrow(model@gebv)))
  }
  
  # PEV
  pev <- compute_pev(model)
  
  # Reliability = 1 - PEV/var_a
  reliability <- 1 - pev / var_a
  reliability[reliability < 0] <- 0
  reliability[reliability > 1] <- 1
  
  return(reliability)
}
