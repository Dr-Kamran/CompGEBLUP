#' Build Additive Relationship Matrix
#'
#' Constructs the additive genomic relationship matrix (A matrix) from marker data
#' using fast C++ implementation with RcppArmadillo.
#'
#' @param M Marker matrix (individuals x markers) coded as -1, 0, 1
#' @param method Method for constructing A matrix ("VanRaden", "UAR", or "UARadj")
#' @param min.MAF Minimum minor allele frequency for filtering markers
#' @param impute.method Method for imputing missing values ("mean")
#' @param return.imputed Logical, return imputed marker matrix
#' @param max.memory.gb Maximum memory to use in GB (default 8)
#' @param use.cpp Use C++ implementation for speed (default TRUE)
#' @param verbose Print messages (default TRUE)
#'
#' @return Additive relationship matrix, or list with A matrix and imputed markers
#' @export
#' @examples
#' \dontrun{
#' M <- matrix(sample(c(-1, 0, 1), 1000, replace = TRUE), 100, 10)
#' A <- build_A_matrix(M, method = "VanRaden")
#' }
build_A_matrix <- function(M, method = c("VanRaden", "UAR", "UARadj"), 
                           min.MAF = 0.01,
                           impute.method = c("mean"),
                           return.imputed = FALSE,
                           max.memory.gb = 8,
                           use.cpp = TRUE,
                           verbose = TRUE) {
  

  # Match arguments
  method <- match.arg(method)
  impute.method <- match.arg(impute.method)
  
  # Input validation
  if (!is.matrix(M)) {
    stop("M must be a matrix")
  }
  
  if (nrow(M) < 2) {
    stop("M must have at least 2 individuals (rows)")
  }
  
  if (ncol(M) < 2) {
    stop("M must have at least 2 markers (columns)")
  }
  
  # Memory check (use as.numeric to avoid integer overflow for large n)
  n_ind <- nrow(M)
  estimated_memory_gb <- (as.numeric(n_ind)^2 * 8) / (1024^3)  # 8 bytes per double
  if (estimated_memory_gb > max.memory.gb) {
    stop(sprintf(
      "Estimated memory usage (%.2f GB) exceeds limit (%.2f GB). ",
      estimated_memory_gb, max.memory.gb),
      "Increase max.memory.gb or reduce number of individuals.")
  }
  
  # Store original marker names
  ind_names <- rownames(M)
  marker_names <- colnames(M)
  
  if (is.null(ind_names)) {
    ind_names <- paste0("ID", seq_len(nrow(M)))
  }
  
  if (is.null(marker_names)) {
    marker_names <- paste0("M", seq_len(ncol(M)))
  }
  
  # Impute missing values
  M_imputed <- impute_markers(M, method = impute.method)
  
  # Filter by MAF
  M_filtered <- filter_markers_maf(M_imputed, min.MAF = min.MAF, verbose = verbose)
  
  if (ncol(M_filtered) < 2) {
    stop("Too few markers remaining after MAF filtering")
  }
  
  # Calculate A matrix based on method
  if (method == "VanRaden" && use.cpp) {
    # Use fast C++ implementation
    A <- calc_A_matrix_cpp(M_filtered)
    
  } else if (method == "VanRaden") {
    # R fallback (VanRaden 2008 method)
    # Calculate allele frequencies
    p <- colMeans(M_filtered + 1) / 2  # Convert -1,0,1 to 0,1,2 then get freq
    
    # Center marker matrix
    M_centered <- scale(M_filtered, center = TRUE, scale = FALSE)
    
    # A = ZZ' / (2 * sum(p * (1-p)))
    denom <- 2 * sum(p * (1 - p))
    if (denom < 1e-10) {
      # Check how many markers have variation
      n_variable <- sum(p > 0.01 & p < 0.99)
      if (n_variable < 10) {
        stop("Cannot compute relationship matrix: ", n_variable, 
             " markers have variation. Need at least 10 polymorphic markers.")
      }
      warning("Denominator near zero (", signif(denom, 3), 
              "), some markers may be monomorphic. Results may be unreliable.")
      denom <- max(denom, 1e-6)  # Use small value, not 1
    }
    A <- tcrossprod(M_centered) / denom
    
  } else if (method == "UAR") {
    # Unified additive relationship (Yang et al. 2010)
    # A = ZZ' / m
    M_centered <- scale(M_filtered, center = TRUE, scale = FALSE)
    A <- tcrossprod(M_centered) / ncol(M_filtered)
    
  } else if (method == "UARadj") {
    # UAR with adjustment for LD
    # A = ZZ' / (m * mean(var(markers)))
    M_centered <- scale(M_filtered, center = TRUE, scale = FALSE)
    marker_vars <- apply(M_centered, 2, var)
    mean_var <- mean(marker_vars)
    if (mean_var < 1e-10) mean_var <- 1
    denom <- ncol(M_filtered) * mean_var
    A <- tcrossprod(M_centered) / denom
  }
  
  # Add dimnames
  dimnames(A) <- list(ind_names, ind_names)
  
  # Ensure symmetry
  A <- (A + t(A)) / 2
  
  # Restore dimnames for imputed markers
  rownames(M_imputed) <- ind_names
  colnames(M_imputed) <- marker_names
  
  if (return.imputed) {
    return(list(A = A, M = M_imputed))
  } else {
    return(A)
  }
}

#' Build Dominance Relationship Matrix
#'
#' Constructs the dominance genomic relationship matrix (D matrix) from marker data
#' using Vitezica et al. (2013) method with fast C++ implementation.
#'
#' @param M Marker matrix (individuals x markers) coded as -1, 0, 1
#' @param method Method for constructing D matrix ("Vitezica" or "Su")
#' @param min.MAF Minimum minor allele frequency for filtering markers
#' @param impute.method Method for imputing missing values ("mean")
#' @param return.imputed Logical, return imputed marker matrix
#' @param max.memory.gb Maximum memory to use in GB (default 8)
#' @param use.cpp Use C++ implementation for speed (default TRUE)
#' @param verbose Print messages (default TRUE)
#'
#' @return Dominance relationship matrix, or list with D matrix and imputed markers
#' @export
build_D_matrix <- function(M, method = c("Vitezica", "Su"), 
                           min.MAF = 0.01,
                           impute.method = c("mean"),
                           return.imputed = FALSE,
                           max.memory.gb = 8,
                           use.cpp = TRUE,
                           verbose = TRUE) {
  
  # Match arguments
  method <- match.arg(method)
  impute.method <- match.arg(impute.method)
  
  # Input validation
  if (!is.matrix(M)) {
    stop("M must be a matrix")
  }
  
  if (nrow(M) < 2) {
    stop("M must have at least 2 individuals (rows)")
  }
  
  # Memory check (use as.numeric to avoid integer overflow for large n)
  n_ind <- nrow(M)
  estimated_memory_gb <- (as.numeric(n_ind)^2 * 8) / (1024^3)
  if (estimated_memory_gb > max.memory.gb) {
    stop(sprintf(
      "Estimated memory usage (%.2f GB) exceeds limit (%.2f GB). ",
      estimated_memory_gb, max.memory.gb),
      "Increase max.memory.gb or reduce number of individuals.")
  }
  
  # Store names
  ind_names <- rownames(M)
  marker_names <- colnames(M)
  
  if (is.null(ind_names)) {
    ind_names <- paste0("ID", seq_len(nrow(M)))
  }
  
  if (is.null(marker_names)) {
    marker_names <- paste0("M", seq_len(ncol(M)))
  }
  
  # Impute missing values
  M_imputed <- impute_markers(M, method = impute.method)
  
  # Filter by MAF
  M_filtered <- filter_markers_maf(M_imputed, min.MAF = min.MAF, verbose = verbose)
  
  if (ncol(M_filtered) < 2) {
    stop("Too few markers remaining after MAF filtering")
  }
  
  if (method == "Vitezica" && use.cpp) {
    # Use fast C++ implementation
    D <- calc_D_matrix_cpp(M_filtered)
    
  } else if (method == "Vitezica") {
    # R fallback - Vitezica et al. (2013) method
    # Calculate allele frequencies
    p <- colMeans(M_filtered + 1) / 2
    
    # Construct dominance marker matrix
    # W = 1 if heterozygote (0), 0 otherwise (-1 or 1)
    W <- matrix(0, nrow = nrow(M_filtered), ncol = ncol(M_filtered))
    W[M_filtered == 0] <- 1
    
    # Center: W - 2*p*(1-p)
    expected_het <- 2 * p * (1 - p)
    W_centered <- sweep(W, 2, expected_het, "-")
    
    # Scale by sum of expected variance
    denom <- sum(expected_het^2)
    if (denom < 1e-10) denom <- 1
    D <- tcrossprod(W_centered) / denom
    
  } else if (method == "Su") {
    # Su et al. (2012) method
    p <- colMeans(M_filtered + 1) / 2
    q <- 1 - p
    
    W <- matrix(0, nrow = nrow(M_filtered), ncol = ncol(M_filtered))
    W[M_filtered == 0] <- 1
    
    expected_het <- 2 * p * q
    W_centered <- sweep(W, 2, expected_het, "-")
    
    # Different scaling
    denom <- 2 * sum(p * q * (1 - 2 * p * q))
    if (denom < 1e-10) denom <- 1
    D <- tcrossprod(W_centered) / denom
  }
  
  # Add dimnames
  dimnames(D) <- list(ind_names, ind_names)
  
  # Ensure symmetry
  D <- (D + t(D)) / 2
  
  # Restore dimnames for imputed markers
  rownames(M_imputed) <- ind_names
  colnames(M_imputed) <- marker_names
  
  if (return.imputed) {
    return(list(D = D, M = M_imputed))
  } else {
    return(D)
  }
}

#' Build Both Additive and Dominance Matrices in One Pass
#'
#' Efficiently computes both A and D relationship matrices simultaneously,
#' avoiding redundant MAF filtering and marker centering.
#'
#' @param M Marker matrix (individuals x markers) coded as -1, 0, 1
#' @param min.MAF Minimum minor allele frequency for filtering markers
#' @param impute.method Method for imputing missing values ("mean")
#' @param use.cpp Use C++ implementation for speed (default TRUE)
#' @param verbose Print messages (default TRUE)
#'
#' @return List with A (additive) and D (dominance) matrices
#' @export
#' @examples
#' \dontrun{
#' M <- matrix(sample(c(-1, 0, 1), 1000, replace = TRUE), 100, 10)
#' matrices <- build_AD_matrices(M)
#' A <- matrices$A
#' D <- matrices$D
#' }
build_AD_matrices <- function(M, 
                               min.MAF = 0.01,
                               impute.method = c("mean"),
                               use.cpp = TRUE,
                               verbose = TRUE) {
  
  impute.method <- match.arg(impute.method)
  
  # Input validation
  if (!is.matrix(M)) {
    stop("M must be a matrix")
  }
  
  # Store original names
  ind_names <- rownames(M)
  marker_names <- colnames(M)
  
  n_ind <- nrow(M)
  n_snp <- ncol(M)
  
  if (verbose) {
    message("Computing A and D matrices: ", n_ind, " individuals x ", n_snp, " markers")
  }
  
  # Impute missing values ONCE
  if (any(is.na(M))) {
    col_means <- colMeans(M, na.rm = TRUE)
    # Handle all-NA columns by setting to 0 (will be filtered by MAF anyway)
    col_means[is.na(col_means)] <- 0
    for (j in seq_len(n_snp)) {
      na_idx <- is.na(M[, j])
      if (any(na_idx)) {
        M[na_idx, j] <- col_means[j]
      }
    }
  }
  
  # Calculate allele frequencies ONCE
  p <- (colMeans(M) + 1) / 2
  q <- 1 - p
  maf <- pmin(p, q)
  
  # Filter by MAF ONCE
  keep <- maf >= min.MAF
  if (sum(keep) < 10) {
    warning("Fewer than 10 markers pass MAF filter. Using all markers.")
    keep <- rep(TRUE, n_snp)
  }
  
  M_filtered <- M[, keep, drop = FALSE]
  p <- p[keep]
  q <- q[keep]
  
  if (verbose && sum(!keep) > 0) {
    message("Filtered ", sum(!keep), " markers with MAF < ", min.MAF)
  }
  
  # ==== ADDITIVE MATRIX (VanRaden) ====
  # Center by 2p
  M_centered <- sweep(M_filtered, 2, 2 * p, "-")
  
  # VanRaden denominator
  denom_A <- 2 * sum(p * q)
  if (denom_A < 1e-6) {
    warning("VanRaden denominator near zero. Using fallback.")
    denom_A <- max(denom_A, 1e-6)
  }
  
  A <- tcrossprod(M_centered) / denom_A
  
  # ==== DOMINANCE MATRIX (Vitezica) ====
  # Dominance coding: heterozygotes = 1, homozygotes = 0
  W <- matrix(0, nrow = n_ind, ncol = sum(keep))
  W[M_filtered == 0] <- 1
  
  expected_het <- 2 * p * q
  W_centered <- sweep(W, 2, expected_het, "-")
  
  denom_D <- sum((2 * p * q)^2)
  if (denom_D < 1e-6) {
    warning("Dominance denominator near zero. Using fallback.")
    denom_D <- max(denom_D, 1e-6)
  }
  
  D <- tcrossprod(W_centered) / denom_D
  
  # Add dimnames
  dimnames(A) <- list(ind_names, ind_names)
  dimnames(D) <- list(ind_names, ind_names)
  
  # Ensure symmetry
  A <- (A + t(A)) / 2
  D <- (D + t(D)) / 2
  
  if (verbose) {
    message("A matrix: diagonal mean = ", round(mean(diag(A)), 3))
    message("D matrix: diagonal mean = ", round(mean(diag(D)), 3))
  }
  
  return(list(A = A, D = D))
}

#' Build Epistatic Relationship Matrix
#'
#' Constructs epistatic relationship matrices (A#A, A#D, or D#D) using
#' Hadamard (element-wise) product with fast C++ implementation.
#'
#' @param M Marker matrix (individuals x markers) coded as -1, 0, 1
#' @param type Type of epistasis ("A#A", "A#D", or "D#D")
#' @param A Optional pre-computed A matrix
#' @param D Optional pre-computed D matrix
#' @param min.MAF Minimum minor allele frequency for filtering markers
#' @param impute.method Method for imputing missing values ("mean")
#' @param use.cpp Use C++ implementation for speed (default TRUE)
#'
#' @return Epistatic relationship matrix
#' @export
build_E_matrix <- function(M, type = c("A#A", "A#D", "D#D"), 
                           A = NULL, D = NULL,
                           min.MAF = 0.01, 
                           impute.method = c("mean"),
                           use.cpp = TRUE) {
  
  # Match arguments
  type <- match.arg(type)
  impute.method <- match.arg(impute.method)
  
  ind_names <- rownames(M)
  if (is.null(ind_names)) {
    ind_names <- paste0("ID", seq_len(nrow(M)))
  }
  
  # Build required matrices if not provided
  if (type %in% c("A#A", "A#D") && is.null(A)) {
    A <- build_A_matrix(M, min.MAF = min.MAF, impute.method = impute.method,
                        use.cpp = use.cpp)
  }
  
  if (type %in% c("A#D", "D#D") && is.null(D)) {
    D <- build_D_matrix(M, min.MAF = min.MAF, impute.method = impute.method,
                        use.cpp = use.cpp)
  }
  
  # Compute Hadamard product
  if (use.cpp) {
    if (type == "A#A") {
      E <- calc_epistatic_matrix_cpp(A, A)
    } else if (type == "A#D") {
      E <- calc_epistatic_matrix_cpp(A, D)
    } else if (type == "D#D") {
      E <- calc_epistatic_matrix_cpp(D, D)
    }
  } else {
    # R fallback
    if (type == "A#A") {
      E <- A * A  # Element-wise multiplication
    } else if (type == "A#D") {
      E <- A * D
    } else if (type == "D#D") {
      E <- D * D
    }
  }
  
  # Add dimnames
  dimnames(E) <- list(ind_names, ind_names)
  
  # Ensure symmetry
  E <- (E + t(E)) / 2
  
  return(E)
}

#' Build G×E Relationship Matrix
#'
#' Constructs relationship matrices for G×E interactions using Kronecker product.
#'
#' @param K Base relationship matrix (A or D)
#' @param n_env Number of environments
#' @param env_cor Optional correlation matrix between environments
#' @param env_names Optional character vector of environment names
#'
#' @return G×E relationship matrix (Kronecker product)
#' @export
build_GE_matrix <- function(K, n_env, env_cor = NULL, env_names = NULL) {
  
  if (!is.matrix(K)) {
    stop("K must be a matrix")
  }
  
  if (nrow(K) != ncol(K)) {
    stop("K must be a square matrix")
  }
  
  if (n_env < 2) {
    stop("n_env must be at least 2")
  }
  
  # If no correlation structure provided, use identity
  if (is.null(env_cor)) {
    env_cor <- diag(n_env)
  }
  
  # Validate env_cor
  if (!is.matrix(env_cor) || nrow(env_cor) != n_env || ncol(env_cor) != n_env) {
    stop("env_cor must be a square matrix of dimension n_env x n_env")
  }
  
  # Validate env_cor/env_names consistency if both provided
  if (!is.null(env_names) && !is.null(rownames(env_cor))) {
    if (!all(rownames(env_cor) == env_names)) {
      warning("env_names does not match rownames(env_cor). Using env_names.")
    }
  }
  
  # Memory check for Kronecker product (use as.numeric to avoid integer overflow)
  n_ind <- nrow(K)
  result_dim <- as.numeric(n_ind) * as.numeric(n_env)
  estimated_memory_gb <- (result_dim^2 * 8) / (1024^3)
  
  if (estimated_memory_gb > 32) {
    stop("Kronecker product would create a ", round(estimated_memory_gb, 2),
         " GB matrix. This exceeds safe memory limits.", call. = FALSE)
  }
  if (estimated_memory_gb > 8) {  # Default 8GB limit
    warning("Kronecker product will create a ", round(estimated_memory_gb, 2), 
            " GB matrix (", format(result_dim, big.mark = ","), " x ", 
            format(result_dim, big.mark = ","), "). ",
            "Consider reducing n_env or n_ind.", call. = FALSE)
  }
  
  # Kronecker product: K ⊗ E
  K_GE <- kronecker(env_cor, K)
  
  # Create dimnames
  ind_names <- rownames(K)
  if (is.null(ind_names)) {
    ind_names <- paste0("ID", seq_len(nrow(K)))
  }
  
  if (is.null(env_names)) {
    env_names <- rownames(env_cor)
  }
  if (is.null(env_names)) {
    env_names <- paste0("E", seq_len(n_env))
  }
  
  # Combine names: ID.Env (consistent with phenotype data format)
  combined_names <- as.vector(outer(ind_names, env_names, paste, sep = "."))
  dimnames(K_GE) <- list(combined_names, combined_names)
  
  return(K_GE)
}

#' Impute Missing Marker Values
#'
#' Imputes missing values in marker matrix using specified method.
#'
#' @param M Marker matrix with possible missing values
#' @param method Imputation method ("mean")
#'
#' @return Imputed marker matrix
#' @keywords internal
impute_markers <- function(M, method = c("mean")) {
  
  method <- match.arg(method)
  
  if (!any(is.na(M))) {
    return(M)
  }
  
  M_imputed <- M
  
  if (method == "mean") {
    # Vectorized mean imputation per marker
    col_means <- colMeans(M, na.rm = TRUE)
    col_means[is.na(col_means)] <- 0  # If all NA, use 0
    
    for (j in seq_len(ncol(M))) {
      missing <- is.na(M[, j])
      if (any(missing)) {
        M_imputed[missing, j] <- col_means[j]
      }
    }
  }
  
  return(M_imputed)
}

#' Filter Markers by Minor Allele Frequency
#'
#' @param M Marker matrix
#' @param min.MAF Minimum minor allele frequency
#'
#' @return Filtered marker matrix
#' @keywords internal
filter_markers_maf <- function(M, min.MAF = 0.01, verbose = TRUE) {
  
  if (min.MAF <= 0) {
    return(M)
  }
  
  # Calculate allele frequencies
  # M is coded as -1, 0, 1
  # Convert to 0, 1, 2 for frequency calculation
  p <- colMeans(M + 1) / 2
  
  # Calculate MAF
  maf <- pmin(p, 1 - p)
  
  # Filter markers
  keep <- maf >= min.MAF
  
  if (sum(keep) == 0) {
    stop("No markers pass MAF filter with threshold ", min.MAF)
  }
  
  M_filtered <- M[, keep, drop = FALSE]
  
  n_removed <- ncol(M) - ncol(M_filtered)
  if (n_removed > 0 && verbose) {
    message("Removed ", n_removed, " markers with MAF < ", min.MAF)
  }
  
  return(M_filtered)
}

#' Make Matrix Positive Definite
#'
#' Ensures a matrix is positive definite. Two modes of operation:
#' \itemize{
#'   \item \code{tol}: Floors eigenvalues below \code{tol} (eigenvalue reconstruction).
#'     This changes the structure of the matrix.
#'   \item \code{ridge}: Adds \code{ridge * I} to the diagonal (diagonal ridge).
#'     This preserves off-diagonal structure and is the standard regularisation approach.
#' }
#' If both are provided, diagonal ridge is applied first, then eigenvalue flooring
#' as a safety net.
#'
#' @param K Square matrix
#' @param tol Tolerance for smallest eigenvalue (eigenvalue floor approach)
#' @param ridge Ridge parameter to add to diagonal. Adds \code{ridge * I} to the
#'   matrix. This is DIFFERENT from eigenvalue flooring via \code{tol}.
#' @param use.cpp Use C++ implementation (default TRUE, only used for eigenvalue floor)
#'
#' @return Positive definite matrix
#' @export
make_positive_definite <- function(K, tol = 1e-6, ridge = NULL, use.cpp = TRUE) {
  
  if (!is.matrix(K)) {
    stop("K must be a matrix")
  }
  
  if (nrow(K) != ncol(K)) {
    stop("K must be a square matrix")
  }
  
  # BUG FIX (v2.3.1): When ridge is provided, add diagonal ridge (K + ridge*I)
  # rather than treating it as eigenvalue floor. Previous versions silently
  # converted ridge to tol, which is a fundamentally different operation.
  if (!is.null(ridge)) {
    K <- K + diag(ridge, nrow(K))
    # After diagonal ridge, check if PD. If so, return early.
    min_eig <- min(eigen(K, symmetric = TRUE, only.values = TRUE)$values)
    if (min_eig >= tol) {
      dimnames_orig <- dimnames(K)
      K <- (K + t(K)) / 2  # Ensure symmetry
      dimnames(K) <- dimnames_orig
      return(K)
    }
    # If still not PD after ridge, fall through to eigenvalue floor
    message("Matrix not PD after diagonal ridge (", ridge,
            "); applying eigenvalue floor (", tol, ") as safety net.")
  }
  
  if (use.cpp) {
    # Use fast C++ implementation
    K_pd <- make_pd_cpp(K, tol)
  } else {
    # R fallback
    # Check if already positive definite
    eig <- eigen(K, symmetric = TRUE)
    min_eig <- min(eig$values)
    
    if (min_eig < tol) {
      # Floor eigenvalues at tol
      eig$values[eig$values < tol] <- tol
      K_pd <- eig$vectors %*% diag(eig$values) %*% t(eig$vectors)
      K_pd <- (K_pd + t(K_pd)) / 2
      message("Adjusted eigenvalues to ensure positive definiteness")
    } else {
      K_pd <- K
    }
  }
  
  # Preserve dimnames
  dimnames(K_pd) <- dimnames(K)
  
  return(K_pd)
}

#' Blend Pedigree and Genomic Relationship Matrices
#'
#' Creates H matrix combining pedigree (A) and genomic (G) information
#'
#' @param A Pedigree-based relationship matrix
#' @param G Genomic relationship matrix
#' @param tau Scaling parameter for G (default: 1)
#' @param omega Scaling parameter for A (default: 1)
#' @param w Weight for G vs A (0-1, default: 0.5)
#'
#' @return H matrix
#' @export
build_H_matrix <- function(A, G, tau = 1, omega = 1, w = 0.5) {
  
  # Input validation
  if (!is.matrix(A) || !is.matrix(G)) {
    stop("A and G must be matrices")
  }
  
  if (nrow(A) != ncol(A) || nrow(G) != ncol(G)) {
    stop("A and G must be square matrices")
  }
  
  stopifnot("w must be between 0 and 1" = w >= 0 && w <= 1)
  
  # Get genotyped individuals (in G)
  gen_ids <- rownames(G)
  all_ids <- rownames(A)
  
  if (is.null(gen_ids) || is.null(all_ids)) {
    stop("A and G must have rownames")
  }
  
  # Check that genotyped individuals are in pedigree
  if (!all(gen_ids %in% all_ids)) {
    stop("All genotyped individuals must be in pedigree matrix A")
  }
  
  # Scale G
  G_scaled <- tau * G
  
  # Extract A_22 (genotyped individuals)
  A_22 <- A[gen_ids, gen_ids]
  
  # Calculate G - A_22
  G_diff <- G_scaled - omega * A_22
  
  # Initialize H with A
  H <- A
  
  # Update H for genotyped individuals
  H[gen_ids, gen_ids] <- H[gen_ids, gen_ids] + w * G_diff
  
  # Ensure symmetry
  H <- (H + t(H)) / 2
  
  return(H)
}

#' Calculate Inbreeding Coefficients from Relationship Matrix
#'
#' @param K Relationship matrix
#'
#' @return Vector of inbreeding coefficients
#' @export
calc_inbreeding <- function(K) {
  
  if (!is.matrix(K)) {
    stop("K must be a matrix")
  }
  
  if (nrow(K) != ncol(K)) {
    stop("K must be a square matrix")
  }
  
  # Inbreeding coefficient = (K_ii - 1) for additive matrix
  f <- diag(K) - 1
  
  # Add names
  names(f) <- rownames(K)
  
  return(f)
}

#' Check if Matrix is Positive Definite
#'
#' @param K Square matrix
#' @param tol Tolerance for smallest eigenvalue
#' @param use.cpp Use C++ implementation (default TRUE)
#'
#' @return Logical, TRUE if positive definite
#' @export
is_positive_definite <- function(K, tol = 1e-8, use.cpp = TRUE) {
  
  if (!is.matrix(K) || nrow(K) != ncol(K)) {
    return(FALSE)
  }
  
  if (use.cpp) {
    return(is_positive_definite_cpp(K, tol))
  } else {
    eig <- eigen(K, symmetric = TRUE, only.values = TRUE)
    return(min(eig$values) > tol)
  }
}
