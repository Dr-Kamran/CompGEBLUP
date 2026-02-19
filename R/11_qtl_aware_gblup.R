#' QTL-Aware Genomic Prediction
#'
#' Functions for genomic prediction that separate major QTL effects (fixed)
#' from polygenic background effects (random). This approach can improve
#' prediction accuracy when major genes are known.
#'
#' @name qtl_aware_gblup
#' @references
#' VanRaden PM (2008). Efficient methods to compute genomic predictions. Journal of Dairy Science 91:4414-4423.
NULL

#' Build QTL-Aware Matrices
#'
#' Constructs relationship matrices that exclude specified QTL markers,
#' and extracts QTL genotypes for use as fixed effects.
#'
#' @param M Marker matrix (individuals x markers) coded as -1, 0, 1
#' @param qtl_markers Character vector of QTL marker names to treat as fixed effects
#' @param qtl_env_markers Optional list of environment-specific QTL markers
#' @param env_names Optional character vector of environment names
#' @param method Method for A matrix ("VanRaden", "UAR")
#' @param min.MAF Minimum MAF for filtering non-QTL markers
#'
#' @return List containing:
#'   \item{Ka}{Reduced additive matrix (excludes QTL markers)}
#'   \item{Xa}{QTL marker genotypes for fixed effects}
#'   \item{A_full}{Full additive matrix (all markers)}
#'   \item{qtl_idx}{Indices of QTL markers}
#'   \item{non_qtl_idx}{Indices of non-QTL markers}
#'
#' @export
#' @examples
#' \dontrun{
#' # After GWAS, get significant markers
#' sig_markers <- get_significant_markers(gwas_res)$marker
#' 
#' # Build QTL-aware matrices
#' qtl_mats <- build_qtl_aware_matrices(M, qtl_markers = sig_markers)
#' 
#' # Fit model with QTLs as fixed effects
#' model <- fit_gblup_qtl(gdata, Ka = qtl_mats$Ka, Xa = qtl_mats$Xa)
#' }
build_qtl_aware_matrices <- function(M, qtl_markers = NULL, 
                                      qtl_env_markers = NULL,
                                      env_names = NULL,
                                      method = c("VanRaden", "UAR"),
                                      min.MAF = 0.01) {
  

  method <- match.arg(method)
  
  if (!is.matrix(M)) {
    stop("M must be a matrix")
  }
  
  # Get marker names
  marker_names <- colnames(M)
  if (is.null(marker_names)) {
    marker_names <- paste0("M", seq_len(ncol(M)))
    colnames(M) <- marker_names
  }
  
  ind_names <- rownames(M)
  if (is.null(ind_names)) {
    ind_names <- paste0("ID", seq_len(nrow(M)))
    rownames(M) <- ind_names
  }
  
  # Validate QTL markers
  if (!is.null(qtl_markers)) {
    missing_qtl <- setdiff(qtl_markers, marker_names)
    if (length(missing_qtl) > 0) {
      warning("QTL markers not found in M: ", 
              paste(head(missing_qtl, 5), collapse = ", "),
              if (length(missing_qtl) > 5) paste0(" (and ", length(missing_qtl) - 5, " more)"))
      qtl_markers <- intersect(qtl_markers, marker_names)
    }
  }
  
  # Get QTL indices
  if (!is.null(qtl_markers) && length(qtl_markers) > 0) {
    qtl_idx <- which(marker_names %in% qtl_markers)
    non_qtl_idx <- which(!marker_names %in% qtl_markers)
    has_qtl <- TRUE
  } else {
    qtl_idx <- integer(0)
    non_qtl_idx <- seq_len(ncol(M))
    has_qtl <- FALSE
  }
  
  # Build full A matrix (all markers)
  A_full <- build_A_matrix(M, method = method, min.MAF = 0, 
                           return.imputed = FALSE, use.cpp = TRUE)
  
  # Build reduced A matrix (excludes QTL markers)
  if (has_qtl && length(non_qtl_idx) > 10) {
    M_reduced <- M[, non_qtl_idx, drop = FALSE]
    Ka <- build_A_matrix(M_reduced, method = method, min.MAF = min.MAF,
                         return.imputed = FALSE, use.cpp = TRUE)
  } else {
    Ka <- A_full
    if (has_qtl) {
      message("Too few non-QTL markers, using full A matrix")
    }
  }
  
  # Extract QTL genotypes for fixed effects
  if (has_qtl) {
    Xa <- M[, qtl_idx, drop = FALSE]
    # Impute missing values in Xa
    for (j in seq_len(ncol(Xa))) {
      na_idx <- is.na(Xa[, j])
      if (any(na_idx)) {
        Xa[na_idx, j] <- mean(Xa[, j], na.rm = TRUE)
      }
    }
    colnames(Xa) <- paste0(marker_names[qtl_idx], "_A")
  } else {
    Xa <- NULL
  }
  
  # Build environment-specific QTL matrices if requested
  KaeE1 <- NULL
  KaeE2 <- NULL
  
  if (!is.null(qtl_env_markers) && !is.null(env_names)) {
    env_result <- build_env_qtl_matrices(M, qtl_env_markers, env_names, Ka)
    KaeE1 <- env_result$KaeE1
    KaeE2 <- env_result$KaeE2
  }
  
  return(list(
    Ka = Ka,
    Xa = Xa,
    A_full = A_full,
    KaeE1 = KaeE1,
    KaeE2 = KaeE2,
    qtl_markers = qtl_markers,
    qtl_idx = qtl_idx,
    non_qtl_idx = non_qtl_idx,
    n_qtl = length(qtl_idx),
    n_markers = ncol(M),
    n_non_qtl = length(non_qtl_idx)
  ))
}

#' Build Environment-Specific QTL Matrices
#'
#' Creates separate kinship matrices for environment-specific QTL effects.
#'
#' @param M Marker matrix
#' @param qtl_env_markers List or data.frame with columns: marker, env
#' @param env_names Character vector of environment names
#' @param Ka Reduced additive matrix (optional, computed if NULL)
#'
#' @return List with KaeE1 (list of per-QTL GxE matrices) and KaeE2 (residual GxE)
#' @keywords internal
build_env_qtl_matrices <- function(M, qtl_env_markers, env_names, Ka = NULL) {
  
  n_ind <- nrow(M)
  n_env <- length(env_names)
  ind_names <- rownames(M)
  marker_names <- colnames(M)
  
  # Parse qtl_env_markers
  if (is.data.frame(qtl_env_markers)) {
    qtl_env_df <- qtl_env_markers
  } else if (is.list(qtl_env_markers)) {
    # Convert list to data frame
    qtl_env_df <- data.frame(
      marker = unlist(qtl_env_markers),
      stringsAsFactors = FALSE
    )
  } else {
    stop("qtl_env_markers must be a data.frame or list")
  }
  
  # Get unique env-specific QTL markers
  env_qtl_markers <- unique(qtl_env_df$marker)
  env_qtl_markers <- intersect(env_qtl_markers, marker_names)
  
  if (length(env_qtl_markers) == 0) {
    return(list(KaeE1 = NULL, KaeE2 = NULL))
  }
  
  # Environment covariance matrix (identity)
  E <- diag(n_env)
  rownames(E) <- colnames(E) <- env_names
  
  # KaeE1: List of per-QTL GxE matrices
  KaeE1 <- vector("list", length(env_qtl_markers))
  names(KaeE1) <- env_qtl_markers
  
  for (i in seq_along(env_qtl_markers)) {
    marker <- env_qtl_markers[i]
    marker_geno <- M[, marker]
    
    # Impute missing
    marker_geno[is.na(marker_geno)] <- mean(marker_geno, na.rm = TRUE)
    
    # Single-marker kinship: xx'
    K_marker <- tcrossprod(marker_geno)
    rownames(K_marker) <- colnames(K_marker) <- ind_names
    
    # Kronecker with environment
    KaeE1[[i]] <- kronecker(E, K_marker, make.dimnames = TRUE)
  }
  
  # KaeE2: Residual GxE matrix (excludes env-specific QTLs)
  non_env_qtl_idx <- which(!marker_names %in% env_qtl_markers)
  
  if (length(non_env_qtl_idx) > 10) {
    M_residual <- M[, non_env_qtl_idx, drop = FALSE]
    Kae_residual <- build_A_matrix(M_residual, min.MAF = 0.01, use.cpp = TRUE)
  } else if (!is.null(Ka)) {
    Kae_residual <- Ka
  } else {
    Kae_residual <- build_A_matrix(M, min.MAF = 0.01, use.cpp = TRUE)
  }
  
  KaeE2 <- kronecker(E, Kae_residual, make.dimnames = TRUE)
  
  return(list(
    KaeE1 = KaeE1,
    KaeE2 = KaeE2,
    env_qtl_markers = env_qtl_markers
  ))
}

#' Fit GBLUP Model with QTL as Fixed Effects
#'
#' Fits a genomic prediction model where major QTLs are treated as fixed effects
#' and the polygenic background is modeled with a reduced kinship matrix.
#'
#' @param gdata GBLUPData object
#' @param Ka Reduced additive matrix (from build_qtl_aware_matrices)
#' @param Xa QTL genotype matrix for fixed effects (from build_qtl_aware_matrices)
#' @param effects Character vector of random effects ("A", "ENV", "AE")
#' @param KaeE1 Optional list of per-QTL GxE matrices
#' @param KaeE2 Optional residual GxE matrix
#' @param ridge Ridge parameter for numerical stability
#' @param nIters Maximum iterations for REML (default: 100). QTL-aware models are
#'   more complex and may need more iterations than standard GBLUP.
#' @param tolParConvLL Convergence tolerance (default: 1e-4)
#' @param verbose Print progress
#'
#' @return GBLUPModel object with additional QTL effect estimates
#' @export
#' @examples
#' \dontrun{
#' # Step 1: Run GWAS
#' K_full <- calc_relationship_matrices(gdata, matrices = "A")
#' gwas_res <- gwas(gdata, K_full)
#' 
#' # Step 2: Get significant QTLs
#' sig_qtl <- get_significant_markers(gwas_res, method = "bonferroni")
#' 
#' # Step 3: Build QTL-aware matrices
#' M <- gdata@genotypes
#' qtl_mats <- build_qtl_aware_matrices(M, qtl_markers = sig_qtl$marker)
#' 
#' # Step 4: Fit model with QTLs as fixed
#' model <- fit_gblup_qtl(gdata, Ka = qtl_mats$Ka, Xa = qtl_mats$Xa)
#' 
#' # Access QTL effects
#' print(model@qtl_effects)
#' }
fit_gblup_qtl <- function(gdata, Ka, Xa = NULL, 
                          effects = c("A", "ENV"),
                          KaeE1 = NULL, KaeE2 = NULL,
                          ridge = 0.001, nIters = 100,
                          tolParConvLL = 1e-4, verbose = TRUE) {
 
  if (!inherits(gdata, "GBLUPData")) {
    stop("gdata must be a GBLUPData object")
  }
  
  pheno <- gdata@phenotypes
  pheno$GID <- as.character(pheno$GID)
  pheno$ENV <- as.character(pheno$ENV)
  
  # Get complete cases
  complete_idx <- !is.na(pheno$trait)
  pheno_complete <- pheno[complete_idx, ]
  y <- as.matrix(pheno_complete$trait)
  
  n_obs <- nrow(pheno_complete)
  
  if (verbose) {
    cat("Fitting QTL-aware GBLUP model\n")
    cat("  Observations:", n_obs, "\n")
    if (!is.null(Xa)) {
      cat("  QTL markers (fixed):", ncol(Xa), "\n")
    }
    cat("  Random effects:", paste(effects, collapse = ", "), "\n")
  }
  
  # Build fixed effects matrix
  # X = [1, Xa] where Xa are QTL genotypes
  X_intercept <- matrix(1, nrow = n_obs, ncol = 1)
  colnames(X_intercept) <- "Intercept"
  
  if (!is.null(Xa) && ncol(Xa) > 0) {
    # Align Xa with phenotype data
    Xa_aligned <- Xa[as.character(pheno_complete$GID), , drop = FALSE]
    
    # Check for NA
    if (any(is.na(Xa_aligned))) {
      # Impute with column means
      for (j in seq_len(ncol(Xa_aligned))) {
        na_idx <- is.na(Xa_aligned[, j])
        if (any(na_idx)) {
          Xa_aligned[na_idx, j] <- mean(Xa_aligned[, j], na.rm = TRUE)
        }
      }
    }
    
    X <- cbind(X_intercept, Xa_aligned)
    has_qtl <- TRUE
  } else {
    X <- X_intercept
    has_qtl <- FALSE
  }
  
  # Build K matrices list
  K_matrices <- list()
  
  # Reduced additive matrix
  if ("A" %in% effects) {
    Ka_reg <- Ka + diag(ridge, nrow(Ka))
    K_matrices$A <- Ka_reg
  }
  
  # Environment matrix (identity)
  if ("ENV" %in% effects) {
    env_levels <- sort(unique(pheno_complete$ENV))
    n_env <- length(env_levels)
    K_env <- diag(n_env)
    rownames(K_env) <- colnames(K_env) <- env_levels
    K_matrices$ENV <- K_env
  }
  
  # GxE matrices
  if ("AE" %in% effects) {
    if (!is.null(KaeE2)) {
      # Use provided residual GxE matrix
      K_matrices$AE <- KaeE2 + diag(ridge, nrow(KaeE2))
    } else {
      # Build standard AE matrix
      env_levels <- sort(unique(pheno_complete$ENV))
      n_env <- length(env_levels)
      E <- diag(n_env)
      rownames(E) <- colnames(E) <- env_levels
      AE <- kronecker(E, Ka, make.dimnames = TRUE)
      K_matrices$AE <- AE + diag(ridge, nrow(AE))
    }
  }
  
  # Add environment-specific QTL effects (KaeE1)
  if (!is.null(KaeE1) && length(KaeE1) > 0) {
    for (i in seq_along(KaeE1)) {
      qtl_name <- names(KaeE1)[i]
      K_matrices[[paste0("QTL_", qtl_name)]] <- KaeE1[[i]] + diag(ridge, nrow(KaeE1[[i]]))
    }
    effects <- c(effects, paste0("QTL_", names(KaeE1)))
  }
  
  # Build design matrices
  design <- build_design_matrices_qtl(pheno_complete, K_matrices, effects)
  
  # Fit model
  fit <- fit_mme(
    y = y,
    X = X,
    Z_list = design$Z_list,
    K_list = design$K_list,
    nIters = nIters,
    tolParConvLL = tolParConvLL,
    tolParInv = ridge,
    verbose = verbose,
    use.emma = (length(design$Z_list) == 1)
  )
  
  # Extract variance components
  theta <- fit$theta
  effect_names_full <- c(design$effect_names, "Residual")
  
  varcomp <- data.frame(
    Component = effect_names_full,
    Variance = theta,
    row.names = NULL,
    stringsAsFactors = FALSE
  )
  
  # Extract QTL effects
  qtl_effects <- NULL
  if (has_qtl) {
    beta <- fit$b
    qtl_beta <- beta[-1, , drop = FALSE]  # Exclude intercept
    rownames(qtl_beta) <- colnames(Xa)
    
    qtl_effects <- data.frame(
      marker = gsub("_A$", "", colnames(Xa)),
      effect = as.numeric(qtl_beta),
      stringsAsFactors = FALSE
    )
    
    if (verbose) {
      cat("\nQTL Fixed Effects:\n")
      print(qtl_effects)
    }
  }
  
  # Extract GEBVs
  gebv_df <- extract_gebvs_qtl(fit, pheno, design, X, Ka, has_qtl, qtl_effects, gdata)
  
  if (verbose) {
    cat("\nVariance Components:\n")
    print(varcomp)
    
    # Heritability
    var_a <- theta[which(design$effect_names == "A")]
    var_total <- sum(theta)
    if (length(var_a) > 0) {
      h2 <- var_a / var_total
      cat("\nPolygenic heritability (h2):", round(h2, 3), "\n")
    }
  }
  
  # Create model object
  model <- new("GBLUPModel",
               model = fit,
               data = gdata,
               gebv = gebv_df,
               varcomp = varcomp,
               description = paste("QTL-GBLUP:", paste(effects, collapse = "+"),
                                   if (has_qtl) paste0(" + ", ncol(Xa), " QTL fixed")),
               effects = effects)
  

  # Add QTL effects as attribute
  attr(model, "qtl_effects") <- qtl_effects
  attr(model, "qtl_markers") <- if (has_qtl) colnames(Xa) else NULL
  
  return(model)
}

#' Build Design Matrices for QTL-Aware Model
#' @keywords internal
build_design_matrices_qtl <- function(pheno, K_matrices, effects) {
  
  pheno$GID <- as.character(pheno$GID)
  pheno$ENV <- as.character(pheno$ENV)
  
  n <- nrow(pheno)
  
  Z_list <- list()
  K_list <- list()
  effect_names <- character()
  
  gid_levels <- sort(unique(pheno$GID))
  env_levels <- sort(unique(pheno$ENV))
  
  for (eff in effects) {
    if (eff == "A" && "A" %in% names(K_matrices)) {
      Z_A <- build_Z_matrix(pheno$GID, gid_levels, sparse = FALSE)
      K_A <- align_K_matrix(K_matrices$A, gid_levels)
      Z_list[[length(Z_list) + 1]] <- Z_A
      K_list[[length(K_list) + 1]] <- K_A
      effect_names <- c(effect_names, "A")
      
    } else if (eff == "ENV" && "ENV" %in% names(K_matrices)) {
      Z_ENV <- build_Z_matrix(pheno$ENV, env_levels, sparse = FALSE)
      K_ENV <- K_matrices$ENV[env_levels, env_levels]
      Z_list[[length(Z_list) + 1]] <- Z_ENV
      K_list[[length(K_list) + 1]] <- K_ENV
      effect_names <- c(effect_names, "ENV")
      
    } else if (eff == "AE" && "AE" %in% names(K_matrices)) {
      gid_env <- paste(pheno$GID, pheno$ENV, sep = ".")
      gid_env_levels <- sort(unique(gid_env))
      Z_AE <- build_Z_matrix(gid_env, gid_env_levels, sparse = FALSE)
      K_AE <- align_K_matrix(K_matrices$AE, gid_env_levels)
      Z_list[[length(Z_list) + 1]] <- Z_AE
      K_list[[length(K_list) + 1]] <- K_AE
      effect_names <- c(effect_names, "AE")
      
    } else if (grepl("^QTL_", eff) && eff %in% names(K_matrices)) {
      gid_env <- paste(pheno$GID, pheno$ENV, sep = ".")
      gid_env_levels <- sort(unique(gid_env))
      Z_QTL <- build_Z_matrix(gid_env, gid_env_levels, sparse = FALSE)
      K_QTL <- align_K_matrix(K_matrices[[eff]], gid_env_levels)
      Z_list[[length(Z_list) + 1]] <- Z_QTL
      K_list[[length(K_list) + 1]] <- K_QTL
      effect_names <- c(effect_names, eff)
    }
  }
  
  return(list(
    Z_list = Z_list,
    K_list = K_list,
    effect_names = effect_names
  ))
}

#' Extract GEBVs from QTL-Aware Model
#' @keywords internal
extract_gebvs_qtl <- function(fit, pheno, design, X, Ka, has_qtl, qtl_effects, gdata = NULL) {
  
  # Find additive effect index
  a_idx <- which(design$effect_names == "A")[1]
  
  if (is.na(a_idx)) {
    warning("No additive effect found")
    return(data.frame(GID = pheno$GID, ENV = pheno$ENV, 
                      GEBV = NA, observed = pheno$trait))
  }
  
  # Get polygenic BLUPs
  u_polygenic <- fit$u[[a_idx]]
  gid_levels <- colnames(design$Z_list[[a_idx]])
  
  if (is.matrix(u_polygenic)) {
    u_polygenic <- as.vector(u_polygenic)
  }
  names(u_polygenic) <- gid_levels
  
  # Calculate total GEBV = polygenic + QTL fixed effects
  pheno$GID <- as.character(pheno$GID)
  pheno$ENV <- as.character(pheno$ENV)
  
  gebv_df <- data.frame(
    GID = pheno$GID,
    ENV = pheno$ENV,
    stringsAsFactors = FALSE
  )
  
  # Add polygenic component
  gebv_df$A_polygenic <- u_polygenic[gebv_df$GID]
  
  # Add QTL fixed effects if present
  if (has_qtl && !is.null(qtl_effects)) {
    # QTL contribution for each individual
    M <- if (!is.null(gdata)) gdata@genotypes else NULL
    if (!is.null(M)) {
      qtl_markers <- qtl_effects$marker
      Xa_full <- M[, qtl_markers, drop = FALSE]
      qtl_contribution <- as.vector(Xa_full %*% qtl_effects$effect)
      names(qtl_contribution) <- rownames(M)
      gebv_df$A_qtl <- qtl_contribution[gebv_df$GID]
    } else {
      gebv_df$A_qtl <- 0
    }
  } else {
    gebv_df$A_qtl <- 0
  }
  
  # Total GEBV
  gebv_df$GEBV <- gebv_df$A_polygenic + gebv_df$A_qtl
  
  # Add observed
  gebv_df$observed <- pheno$trait
  
  # Clean up intermediate columns for output
  gebv_out <- gebv_df[, c("GID", "ENV", "GEBV", "observed")]
  
  # Store detailed breakdown as attribute
  attr(gebv_out, "components") <- gebv_df[, c("GID", "ENV", "A_polygenic", "A_qtl", "GEBV")]
  
  return(gebv_out)
}

#' Cross-Validation for QTL-Aware Model
#'
#' Performs cross-validation with QTL markers as fixed effects.
#'
#' @param gdata GBLUPData object
#' @param qtl_markers Character vector of QTL marker names
#' @param qtl_env_markers Optional environment-specific QTL markers
#' @param scheme CV scheme: "CV1" (new genotypes) or "CV2" (new environments)
#' @param n_folds Number of folds for CV1
#' @param n_reps Number of repetitions
#' @param effects Random effects to include
#' @param verbose Print progress
#'
#' @return CVResult object
#' @export
cv_gblup_qtl <- function(gdata, qtl_markers = NULL,
                         qtl_env_markers = NULL,
                         scheme = c("CV1", "CV2"),
                         n_folds = 5, n_reps = 1,
                         effects = c("A", "ENV"),
                         verbose = TRUE) {
  
  scheme <- match.arg(scheme)
  
  M <- gdata@genotypes
  pheno <- gdata@phenotypes
  
  # Build QTL-aware matrices
  qtl_mats <- build_qtl_aware_matrices(M, qtl_markers = qtl_markers)
  Ka <- qtl_mats$Ka
  Xa <- qtl_mats$Xa
  
  # Get unique IDs/environments
  ids <- unique(pheno$GID)
  envs <- unique(pheno$ENV)
  
  results <- data.frame()
  
  for (rep in seq_len(n_reps)) {
    if (verbose) cat("Repetition", rep, "of", n_reps, "\n")
    
    if (scheme == "CV1") {
      # Leave genotypes out
      set.seed(rep * 123)
      folds <- sample(rep(1:n_folds, length.out = length(ids)))
      names(folds) <- ids
      
      for (fold in seq_len(n_folds)) {
        if (verbose) cat("  Fold", fold, "...\n")
        
        test_ids <- names(folds[folds == fold])
        
        # Mask test phenotypes
        pheno_masked <- pheno
        pheno_masked$trait[pheno_masked$GID %in% test_ids] <- NA
        
        gdata_train <- new("GBLUPData",
                          genotypes = M,
                          phenotypes = pheno_masked,
                          maf = gdata@maf)
        
        # Fit model
        model <- tryCatch({
          fit_gblup_qtl(gdata_train, Ka = Ka, Xa = Xa, 
                        effects = effects, verbose = FALSE)
        }, error = function(e) NULL)
        
        if (!is.null(model)) {
          # Get predictions for test set
          gebv <- model@gebv
          test_obs <- pheno[pheno$GID %in% test_ids, ]
          
          pred_vals <- gebv$GEBV[match(paste(test_obs$GID, test_obs$ENV),
                                       paste(gebv$GID, gebv$ENV))]
          
          acc <- cor(test_obs$trait, pred_vals, use = "complete.obs")
          
          results <- rbind(results, data.frame(
            rep = rep, fold = fold, scheme = scheme,
            predictive_ability = acc, stringsAsFactors = FALSE
          ))
        }
      }
      
    } else if (scheme == "CV2") {
      # Leave environments out
      for (env in envs) {
        if (verbose) cat("  Environment:", env, "\n")
        
        # Mask test environment
        pheno_masked <- pheno
        pheno_masked$trait[pheno_masked$ENV == env] <- NA
        
        gdata_train <- new("GBLUPData",
                          genotypes = M,
                          phenotypes = pheno_masked,
                          maf = gdata@maf)
        
        # Fit model
        model <- tryCatch({
          fit_gblup_qtl(gdata_train, Ka = Ka, Xa = Xa,
                        effects = effects, verbose = FALSE)
        }, error = function(e) NULL)
        
        if (!is.null(model)) {
          gebv <- model@gebv
          test_obs <- pheno[pheno$ENV == env, ]
          
          pred_vals <- gebv$GEBV[match(paste(test_obs$GID, test_obs$ENV),
                                       paste(gebv$GID, gebv$ENV))]
          
          acc <- cor(test_obs$trait, pred_vals, use = "complete.obs")
          
          results <- rbind(results, data.frame(
            rep = rep, fold = env, scheme = scheme,
            predictive_ability = acc, stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # Create CVResult object
  cv_result <- new("CVResult",
                   metrics = results,
                   predictions = data.frame(),
                   scheme = scheme,
                   n_folds = n_folds,
                   effects = effects)
  
  if (verbose) {
    cat("\nCross-validation results:\n")
    cat("  Mean predictive ability:", round(mean(results$predictive_ability, na.rm = TRUE), 3), "\n")
    cat("  SD predictive ability:", round(sd(results$predictive_ability, na.rm = TRUE), 3), "\n")
  }
  
  return(cv_result)
}

#' Automatic QTL Detection and Model Fitting Pipeline
#'
#' Combines GWAS, QTL detection, and QTL-aware model fitting in one function.
#'
#' @param gdata GBLUPData object
#' @param gwas_threshold P-value threshold for QTL detection (default: 1e-4)
#' @param gwas_method Multiple testing correction ("bonferroni", "fdr", "none")
#' @param max_qtl Maximum number of QTLs to include
#' @param effects Random effects to include
#' @param verbose Print progress
#'
#' @return List with GWAS results, QTL markers, and fitted model
#' @export
#' @examples
#' \dontrun{
#' result <- fit_gblup_auto_qtl(gdata, gwas_threshold = 0.05, 
#'                               gwas_method = "bonferroni")
#' print(result$qtl_markers)
#' print(result$model@varcomp)
#' }
fit_gblup_auto_qtl <- function(gdata, gwas_threshold = 1e-4,
                                gwas_method = c("bonferroni", "fdr", "none"),
                                max_qtl = 20,
                                effects = c("A", "ENV"),
                                verbose = TRUE) {
  
  gwas_method <- match.arg(gwas_method)
  
  if (verbose) cat("Step 1: Running GWAS for QTL detection...\n")
  
  # Build full A matrix for GWAS
  K_full <- calc_relationship_matrices(gdata, matrices = "A")
  
  # Run GWAS
  gwas_res <- gwas(gdata, K_matrices = K_full, verbose = FALSE)
  
  # Get significant markers
  if (gwas_method == "bonferroni") {
    threshold <- gwas_threshold / gwas_res@n_markers
    sig_markers <- gwas_res@results[gwas_res@results$p_value < threshold, ]
  } else if (gwas_method == "fdr") {
    sig_markers <- get_significant_markers(gwas_res, method = "fdr", alpha = gwas_threshold)
  } else {
    sig_markers <- gwas_res@results[gwas_res@results$p_value < gwas_threshold, ]
  }
  
  # Limit number of QTLs
  if (nrow(sig_markers) > max_qtl) {
    sig_markers <- sig_markers[order(sig_markers$p_value), ]
    sig_markers <- sig_markers[1:max_qtl, ]
  }
  
  qtl_markers <- sig_markers$marker
  
  if (verbose) {
    cat("  Detected", length(qtl_markers), "significant QTL markers\n")
    if (length(qtl_markers) > 0) {
      cat("  Top QTLs:", paste(head(qtl_markers, 5), collapse = ", "), "\n")
    }
  }
  
  # Build QTL-aware matrices
  if (verbose) cat("\nStep 2: Building QTL-aware matrices...\n")
  
  M <- gdata@genotypes
  qtl_mats <- build_qtl_aware_matrices(M, qtl_markers = qtl_markers)
  
  if (verbose) {
    cat("  Markers in polygenic K:", qtl_mats$n_non_qtl, "\n")
    cat("  QTL markers as fixed:", qtl_mats$n_qtl, "\n")
  }
  
  # Fit model
  if (verbose) cat("\nStep 3: Fitting QTL-aware model...\n")
  
  model <- fit_gblup_qtl(gdata, Ka = qtl_mats$Ka, Xa = qtl_mats$Xa,
                         effects = effects, verbose = verbose)
  
  # Compare with standard model
  if (verbose) {
    cat("\nStep 4: Comparison with standard GBLUP...\n")
    
    model_std <- fit_gblup(gdata, K_full, effects = effects, verbose = FALSE)
    
    acc_qtl <- cor(model@gebv$GEBV, model@gebv$observed, use = "complete.obs")
    acc_std <- cor(model_std@gebv$GEBV, model_std@gebv$observed, use = "complete.obs")
    
    cat("  Standard GBLUP predictive ability:", round(acc_std, 3), "\n")
    cat("  QTL-aware GBLUP predictive ability:", round(acc_qtl, 3), "\n")
    cat("  Improvement:", round((acc_qtl - acc_std) * 100, 1), "%\n")
  }
  
  return(list(
    gwas_result = gwas_res,
    qtl_markers = qtl_markers,
    qtl_effects = attr(model, "qtl_effects"),
    qtl_matrices = qtl_mats,
    model = model
  ))
}

#' Get QTL Effects from Model
#'
#' Extracts the estimated QTL fixed effects from a QTL-aware model.
#'
#' @param model GBLUPModel fitted with fit_gblup_qtl
#' @return Data frame with marker names and effects
#' @export
get_qtl_effects <- function(model) {
  
  qtl_effects <- attr(model, "qtl_effects")
  
  if (is.null(qtl_effects)) {
    message("No QTL effects in model (standard GBLUP?)")
    return(NULL)
  }
  
  return(qtl_effects)
}

#' Summary of GEBV Components
#'
#' Shows the breakdown of GEBV into polygenic and QTL components.
#'
#' @param model GBLUPModel fitted with fit_gblup_qtl
#' @return Data frame with component breakdown
#' @export
get_gebv_components <- function(model) {
  
  components <- attr(model@gebv, "components")
  
  if (is.null(components)) {
    message("No component breakdown available")
    return(model@gebv)
  }
  
  return(components)
}
