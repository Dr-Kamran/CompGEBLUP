#' Design Matrix Construction
#'
#' Functions to build design matrices for mixed models with optional
#' sparse matrix support for memory efficiency.
#'
#' @name design_matrices
#' @keywords internal
NULL

#' Build Design Matrices for GBLUP
#'
#' Constructs X (fixed) and Z (random) design matrices based on
#' the phenotype data and specified effects. Supports sparse matrices
#' for improved memory efficiency with large datasets.
#'
#' Built-in effects: \code{"A"}, \code{"D"}, \code{"ENV"}, \code{"AE"},
#' \code{"DE"}, \code{"AA"}, \code{"AAE"}. Any other effect name is treated
#' as a custom effect — a matching K matrix must be provided in
#' \code{K_matrices}. Custom kernels are auto-detected as GID-level
#' (n_gid x n_gid) or GID.ENV-level (n_gid*n_env x n_gid*n_env) based
#' on dimensions.
#'
#' @param pheno Phenotype data frame with GID, ENV, trait columns
#' @param K_matrices List of relationship matrices. Names must match \code{effects}.
#'   For custom effects, include a matrix with the matching name.
#' @param effects Character vector of effects to include. Can mix built-in
#'   (\code{"A"}, \code{"D"}, etc.) and custom names (e.g. \code{"RN"}, \code{"K_comb"}).
#' @param use.sparse Use sparse matrices for Z (default FALSE, set TRUE for large n)
#' @param fixed Optional fixed effects specification. Can be NULL (intercept only, default),
#'   a one-sided formula (e.g., ~ TESTER), or a pre-built design matrix
#'
#' @return List with X, Z_list, K_list, and effect_names
#' @keywords internal
build_design_matrices <- function(pheno, K_matrices, effects, use.sparse = FALSE, fixed = NULL) {
  
  # Ensure character types
  pheno$GID <- as.character(pheno$GID)
  pheno$ENV <- as.character(pheno$ENV)
  
  n <- nrow(pheno)
  
  # Fixed effects: build X matrix based on input
  if (is.null(fixed)) {
    # Default: intercept only
    X <- matrix(1, nrow = n, ncol = 1)
    colnames(X) <- "(Intercept)"
  } else if (inherits(fixed, "formula")) {
    # Formula provided: use model.matrix
    X <- model.matrix(fixed, data = pheno)
    
    # Check if intercept is missing and add it if needed
    # Look for intercept column by name (handles both "(Intercept)" and "Intercept")
    has_intercept <- any(grepl("^\\(Intercept\\)$|^Intercept$", colnames(X), ignore.case = TRUE))
    if (!has_intercept) {
      X_intercept <- matrix(1, nrow = n, ncol = 1)
      colnames(X_intercept) <- "(Intercept)"
      X <- cbind(X_intercept, X)
    }
  } else if (is.matrix(fixed)) {
    # Pre-built matrix provided
    if (nrow(fixed) != n) {
      stop("Fixed effects matrix must have ", n, " rows to match phenotype data")
    }
    X <- fixed
    
    # Ensure intercept is included by checking column names
    # Look for columns named "Intercept" or "(Intercept)"
    has_intercept <- any(grepl("^\\(Intercept\\)$|^Intercept$", colnames(X), ignore.case = TRUE))
    if (!has_intercept) {
      warning("No intercept column detected in fixed effects matrix. Adding intercept column.",
              call. = FALSE)
      X_intercept <- matrix(1, nrow = n, ncol = 1)
      colnames(X_intercept) <- "(Intercept)"
      X <- cbind(X_intercept, X)
    }
  } else {
    stop("fixed must be NULL, a formula, or a matrix")
  }
  
  # Initialize random effects
  Z_list <- list()
  K_list <- list()
  effect_names <- character()
  
  # Get unique levels
  gid_levels <- sort(unique(pheno$GID))
  env_levels <- sort(unique(pheno$ENV))
  
  # Create environment identity matrix (used for ENV effect and GxE auto-computation)
  n_env <- length(env_levels)
  I_env <- diag(n_env)
  rownames(I_env) <- colnames(I_env) <- env_levels
  
  # Additive genetic effect (A)
  if ("A" %in% effects) {
    if (!("A" %in% names(K_matrices))) {
      stop("Effect 'A' specified but no 'A' matrix in K_matrices")
    }
    
    Z_A <- build_Z_matrix(pheno$GID, gid_levels, sparse = use.sparse)
    K_A <- align_K_matrix(K_matrices$A, gid_levels)
    
    Z_list[[length(Z_list) + 1]] <- Z_A
    K_list[[length(K_list) + 1]] <- K_A
    effect_names <- c(effect_names, "A")
  }
  
  # Dominance effect (D)
  if ("D" %in% effects) {
    if (!("D" %in% names(K_matrices))) {
      stop("Effect 'D' specified but no 'D' matrix in K_matrices")
    }
    
    Z_D <- build_Z_matrix(pheno$GID, gid_levels, sparse = use.sparse)
    K_D <- align_K_matrix(K_matrices$D, gid_levels)
    
    Z_list[[length(Z_list) + 1]] <- Z_D
    K_list[[length(K_list) + 1]] <- K_D
    effect_names <- c(effect_names, "D")
  }
  
  # Environment effect (ENV)
  if ("ENV" %in% effects) {
    Z_ENV <- build_Z_matrix(pheno$ENV, env_levels, sparse = use.sparse)
    
    Z_list[[length(Z_list) + 1]] <- Z_ENV
    K_list[[length(K_list) + 1]] <- I_env
    effect_names <- c(effect_names, "ENV")
  }
  
  # Additive x Environment interaction (AE)
  if ("AE" %in% effects) {
    # Check if AE matrix is available, or can be auto-computed from A
    if (!("AE" %in% names(K_matrices))) {
      if ("A" %in% names(K_matrices)) {
        # Auto-compute AE matrix from A using Kronecker product
        # NOTE (v2.3.1): If K_matrices$A has been pre-regularised (ridge added),
        # the AE matrix will inherit that ridge. For strict separation, pass
        # unridged A in K_matrices and use ridge_per_matrix in fit_gblup to
        # add ridge after design matrix construction.
        message("Auto-computing AE matrix from A matrix using Kronecker product")
        
        # Compute Kronecker product: AE = I_env ⊗ A (environment first, matching build_GE_matrix)
        K_A <- K_matrices$A
        K_AE_computed <- kronecker(I_env, K_A)
        
        # Create dimnames manually to match GID.ENV format
        ind_names <- rownames(K_A)
        if (is.null(ind_names)) ind_names <- paste0("ID", seq_len(nrow(K_A)))
        combined_names <- as.vector(outer(ind_names, env_levels, paste, sep = "."))
        dimnames(K_AE_computed) <- list(combined_names, combined_names)
        
        K_matrices$AE <- K_AE_computed
      } else {
        stop("Effect 'AE' specified but no 'AE' matrix in K_matrices. ",
             "Provide 'AE' in K_matrices, or provide 'A' so it can be auto-computed. ",
             "Use build_GE_matrix() or calc_relationship_matrices() to pre-compute.")
      }
    }
    
    # Create GID:ENV interaction term
    gid_env <- paste(pheno$GID, pheno$ENV, sep = ".")
    gid_env_levels <- sort(unique(gid_env))
    
    Z_AE <- build_Z_matrix(gid_env, gid_env_levels, sparse = use.sparse)
    K_AE <- align_K_matrix(K_matrices$AE, gid_env_levels)
    
    Z_list[[length(Z_list) + 1]] <- Z_AE
    K_list[[length(K_list) + 1]] <- K_AE
    effect_names <- c(effect_names, "AE")
  }
  
  # Dominance x Environment interaction (DE)
  if ("DE" %in% effects) {
    # Check if DE matrix is available, or can be auto-computed from D
    if (!("DE" %in% names(K_matrices))) {
      if ("D" %in% names(K_matrices)) {
        # Auto-compute DE matrix from D using Kronecker product
        message("Auto-computing DE matrix from D matrix using Kronecker product")
        
        # Compute Kronecker product: DE = I_env ⊗ D (environment first, matching build_GE_matrix)
        K_D <- K_matrices$D
        K_DE_computed <- kronecker(I_env, K_D)
        
        # Create dimnames manually to match GID.ENV format
        ind_names <- rownames(K_D)
        if (is.null(ind_names)) ind_names <- paste0("ID", seq_len(nrow(K_D)))
        combined_names <- as.vector(outer(ind_names, env_levels, paste, sep = "."))
        dimnames(K_DE_computed) <- list(combined_names, combined_names)
        
        K_matrices$DE <- K_DE_computed
      } else {
        stop("Effect 'DE' specified but no 'DE' matrix in K_matrices. ",
             "Provide 'DE' in K_matrices, or provide 'D' so it can be auto-computed. ",
             "Use build_GE_matrix() or calc_relationship_matrices() to pre-compute.")
      }
    }
    
    gid_env <- paste(pheno$GID, pheno$ENV, sep = ".")
    gid_env_levels <- sort(unique(gid_env))
    
    Z_DE <- build_Z_matrix(gid_env, gid_env_levels, sparse = use.sparse)
    K_DE <- align_K_matrix(K_matrices$DE, gid_env_levels)
    
    Z_list[[length(Z_list) + 1]] <- Z_DE
    K_list[[length(K_list) + 1]] <- K_DE
    effect_names <- c(effect_names, "DE")
  }
  
  # Epistasis (AA)
  if ("AA" %in% effects) {
    if (!("AA" %in% names(K_matrices))) {
      stop("Effect 'AA' specified but no 'AA' matrix in K_matrices")
    }
    
    Z_AA <- build_Z_matrix(pheno$GID, gid_levels, sparse = use.sparse)
    K_AA <- align_K_matrix(K_matrices$AA, gid_levels)
    
    Z_list[[length(Z_list) + 1]] <- Z_AA
    K_list[[length(K_list) + 1]] <- K_AA
    effect_names <- c(effect_names, "AA")
  }
  
  # Epistasis x Environment (AAE)
  if ("AAE" %in% effects) {
    # Check if AAE matrix is available, or can be auto-computed from AA
    if (!("AAE" %in% names(K_matrices))) {
      if ("AA" %in% names(K_matrices)) {
        # Auto-compute AAE matrix from AA using Kronecker product
        message("Auto-computing AAE matrix from AA matrix using Kronecker product")
        
        # Compute Kronecker product: AAE = I_env ⊗ AA (environment first, matching build_GE_matrix)
        K_AA <- K_matrices$AA
        K_AAE_computed <- kronecker(I_env, K_AA)
        
        # Create dimnames manually to match GID.ENV format
        ind_names <- rownames(K_AA)
        if (is.null(ind_names)) ind_names <- paste0("ID", seq_len(nrow(K_AA)))
        combined_names <- as.vector(outer(ind_names, env_levels, paste, sep = "."))
        dimnames(K_AAE_computed) <- list(combined_names, combined_names)
        
        K_matrices$AAE <- K_AAE_computed
      } else {
        stop("Effect 'AAE' specified but no 'AAE' matrix in K_matrices. ",
             "Provide 'AAE' in K_matrices, or provide 'AA' so it can be auto-computed. ",
             "Use build_GE_matrix() or calc_relationship_matrices() to pre-compute.")
      }
    }
    
    gid_env <- paste(pheno$GID, pheno$ENV, sep = ".")
    gid_env_levels <- sort(unique(gid_env))
    
    Z_AAE <- build_Z_matrix(gid_env, gid_env_levels, sparse = use.sparse)
    K_AAE <- align_K_matrix(K_matrices$AAE, gid_env_levels)
    
    Z_list[[length(Z_list) + 1]] <- Z_AAE
    K_list[[length(K_list) + 1]] <- K_AAE
    effect_names <- c(effect_names, "AAE")
  }
  
  # =========================================================================
  # CUSTOM EFFECTS (v2.3.1): Handle any effect name not in the built-in set.
  # This allows users to pass arbitrary relationship matrices (e.g. reaction
  # norm kernels, combined kernels, factor-analytic kernels) with any name.
  #
  # Convention: If an effect name is not one of the 7 built-in types
  # (A, D, ENV, AE, DE, AA, AAE) but a matching K matrix exists in

  # K_matrices, treat it as a GID-level random effect (same Z as A/D).
  # If the K matrix dimensions match n_gid × n_gid, use GID-level Z.
  # If dimensions match n_gid*n_env × n_gid*n_env, use GID.ENV-level Z
  # (i.e. it's a GxE-type kernel).
  # =========================================================================
  builtin_effects <- c("A", "D", "ENV", "AE", "DE", "AA", "AAE")
  custom_effects <- setdiff(effects, builtin_effects)
  
  for (eff in custom_effects) {
    if (!(eff %in% names(K_matrices))) {
      warning("Custom effect '", eff, "' specified but no matching matrix in K_matrices. ",
              "Skipping this effect.", call. = FALSE)
      next
    }
    
    K_custom <- K_matrices[[eff]]
    n_gid <- length(gid_levels)
    n_gid_env <- n_gid * n_env
    
    if (nrow(K_custom) == n_gid && ncol(K_custom) == n_gid) {
      # GID-level kernel (same structure as A, D, AA)
      Z_custom <- build_Z_matrix(pheno$GID, gid_levels, sparse = use.sparse)
      K_custom <- align_K_matrix(K_custom, gid_levels)
      message("Custom effect '", eff, "': GID-level kernel (", n_gid, " x ", n_gid, ")")
    } else if (nrow(K_custom) == n_gid_env && ncol(K_custom) == n_gid_env) {
      # GID×ENV-level kernel (same structure as AE, DE, AAE)
      gid_env <- paste(pheno$GID, pheno$ENV, sep = ".")
      gid_env_levels <- sort(unique(gid_env))
      Z_custom <- build_Z_matrix(gid_env, gid_env_levels, sparse = use.sparse)
      K_custom <- align_K_matrix(K_custom, gid_env_levels)
      message("Custom effect '", eff, "': GID.ENV-level kernel (",
              n_gid_env, " x ", n_gid_env, ")")
    } else {
      warning("Custom effect '", eff, "': K matrix dimensions (",
              nrow(K_custom), " x ", ncol(K_custom),
              ") do not match n_gid (", n_gid, ") or n_gid*n_env (",
              n_gid_env, "). Skipping.", call. = FALSE)
      next
    }
    
    Z_list[[length(Z_list) + 1]] <- Z_custom
    K_list[[length(K_list) + 1]] <- K_custom
    effect_names <- c(effect_names, eff)
  }
  
  if (length(Z_list) == 0) {
    stop("No random effects specified. Check that effect names match K_matrices names.")
  }
  
  # Convert sparse matrices to dense for MME solver
  # NOTE: Current MME solver requires dense matrices. Sparse storage saves memory
  # during construction but must be converted for solving.
  # TODO: Implement sparse MME solver for large-scale applications (n > 10,000)
  # This would require:
  # 1. Sparse Cholesky factorization for V inverse
  # 2. Sparse-dense matrix products in REML iterations
  # 3. Efficient sparse ZKZ' computation
  Z_list <- lapply(Z_list, function(Z) {
    if (inherits(Z, "sparseMatrix")) {
      # Issue informative message for large matrices
      if (nrow(Z) * ncol(Z) > 1e7) {
        message("Note: Converting large sparse Z matrix to dense. ",
                "For very large datasets (n > 10,000), consider using specialized ",
                "sparse MME software.")
      }
      as.matrix(Z)
    } else {
      Z
    }
  })
  
  return(list(
    X = X,
    Z_list = Z_list,
    K_list = K_list,
    effect_names = effect_names
  ))
}

#' Build Z Matrix (Incidence Matrix)
#'
#' Creates a design matrix mapping observations to levels.
#' Supports both dense and sparse matrix output.
#'
#' @param factor_vector Vector of factor values for each observation
#' @param levels Character vector of all levels (in desired order)
#' @param sparse Logical, return sparse matrix (default FALSE)
#'
#' @return Design matrix (n x q)
#' @keywords internal
build_Z_matrix <- function(factor_vector, levels, sparse = FALSE) {
  
  n <- length(factor_vector)
  q <- length(levels)
  
  if (sparse) {
    # Create sparse matrix using Matrix package
    # Find which level each observation belongs to
    factor_idx <- match(factor_vector, levels)
    
    # Remove NA matches
    valid <- !is.na(factor_idx)
    
    Z <- Matrix::sparseMatrix(
      i = which(valid),
      j = factor_idx[valid],
      x = 1,
      dims = c(n, q),
      dimnames = list(NULL, levels)
    )
  } else {
    # Create dense indicator matrix
    Z <- matrix(0, nrow = n, ncol = q)
    colnames(Z) <- levels
    
    # Vectorized assignment using match
    factor_idx <- match(factor_vector, levels)
    valid <- !is.na(factor_idx)
    
    for (i in which(valid)) {
      Z[i, factor_idx[i]] <- 1
    }
  }
  
  return(Z)
}

#' Build Z Matrix Using C++ (Fast Version)
#'
#' Creates a sparse incidence matrix using C++ for maximum speed.
#'
#' @param factor_vector Vector of factor values for each observation
#' @param levels Character vector of all levels (in desired order)
#'
#' @return Sparse design matrix
#' @keywords internal
build_Z_matrix_cpp <- function(factor_vector, levels) {
  
  factor_idx <- match(factor_vector, levels)
  
  # Handle NAs by setting to 0 (will create empty rows)
  factor_idx[is.na(factor_idx)] <- 0L
  
  # Use C++ function
  Z_sparse <- create_Z_sparse_cpp(as.integer(factor_idx), length(levels))
  
  # Add column names
  colnames(Z_sparse) <- levels
  
  return(Z_sparse)
}

#' Align K Matrix to Factor Levels
#'
#' Ensures the K matrix rows/columns match the specified levels.
#'
#' @param K Relationship matrix
#' @param levels Character vector of required levels
#'
#' @return Aligned K matrix
#' @keywords internal
align_K_matrix <- function(K, levels) {
  
  if (is.null(rownames(K)) || is.null(colnames(K))) {
    if (nrow(K) == length(levels)) {
      rownames(K) <- colnames(K) <- levels
      return(K)
    } else {
      stop("K matrix has no rownames and dimensions don't match levels (",
           nrow(K), " vs ", length(levels), ")")
    }
  }
  
  # Check if all levels are in K
  missing_levels <- setdiff(levels, rownames(K))
  if (length(missing_levels) > 0) {
    stop("Levels missing from K matrix: ", 
         paste(head(missing_levels, 5), collapse = ", "),
         if (length(missing_levels) > 5) paste0(" (and ", length(missing_levels) - 5, " more)"))
  }
  
  # Subset and reorder K to match levels
  K_aligned <- K[levels, levels, drop = FALSE]
  
  return(K_aligned)
}

#' Compute ZKZ' Efficiently
#'
#' Computes Z * K * Z' using either dense or sparse operations
#' depending on input types.
#'
#' @param Z Incidence matrix (can be sparse)
#' @param K Relationship matrix
#'
#' @return ZKZ' matrix
#' @keywords internal
compute_ZKZt <- function(Z, K) {
  
  if (inherits(Z, "sparseMatrix")) {
    # Sparse computation
    ZK <- Z %*% K
    ZKZt <- as.matrix(ZK %*% Matrix::t(Z))
  } else {
    # Dense computation - use C++ if available
    ZKZt <- tryCatch({
      compute_ZKZt_cpp(Z, K)
    }, error = function(e) {
      # R fallback
      Z %*% K %*% t(Z)
    })
  }
  
  return(ZKZt)
}

#' Create Model Matrix from Formula
#'
#' Creates fixed effects design matrix from a formula and data frame.
#'
#' @param formula Formula for fixed effects
#' @param data Data frame
#'
#' @return Design matrix
#' @export
create_model_matrix <- function(formula, data) {
  
  # Use model.matrix with na.pass to handle missing values
  X <- model.matrix(formula, data = data, na.action = na.pass)
  
  return(X)
}

#' Get Factor Level Indices
#'
#' Returns the integer indices mapping observations to factor levels.
#'
#' @param factor_vector Vector of factor values
#' @param levels Optional vector of levels (if NULL, uses unique values)
#'
#' @return Integer vector of level indices (1-based)
#' @keywords internal
get_factor_indices <- function(factor_vector, levels = NULL) {
  
  if (is.null(levels)) {
    levels <- sort(unique(factor_vector))
  }
  
  idx <- match(factor_vector, levels)
  
  return(idx)
}
