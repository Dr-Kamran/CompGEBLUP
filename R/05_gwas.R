#' Genome-Wide Association Study (GWAS)
#'
#' Performs GWAS using EMMAX-style efficient mixed linear model with 
#' spectral decomposition for fast marker testing.
#'
#' @param gdata GBLUPData object
#' @param K_matrices List of relationship matrices for population structure control
#' @param effects Character vector of random effects to include in the model
#' @param min.MAF Minimum minor allele frequency for filtering markers
#' @param P3D Logical, if TRUE variance components estimated once (default TRUE)
#' @param n.PC Number of principal components to include as fixed effects
#' @param method Method for p-value calculation ("Wald" or "LRT")
#' @param n.cores Number of cores for parallel processing (default 1)
#' @param use.cpp Use C++ implementation for speed (default TRUE)
#' @param verbose Logical, print progress messages
#'
#' @return GWASResult object
#' @export
#' @examples
#' \dontrun{
#' gdata <- simulate_genotypes(n_ind = 200, n_snp = 1000)
#' gdata <- simulate_phenotypes(gdata, n_env = 1, h2 = 0.6, n_qtl = 10)
#' K <- calc_relationship_matrices(gdata, matrices = "A")
#' gwas_res <- gwas(gdata, K_matrices = K)
#' plot_manhattan(gwas_res)
#' plot_qq(gwas_res)
#' }
gwas <- function(gdata, K_matrices = NULL, effects = "A",
                 min.MAF = 0.05, P3D = TRUE, n.PC = 0,
                 method = c("Wald", "LRT"), n.cores = 1,
                 use.cpp = TRUE, verbose = TRUE) {
  
  method <- match.arg(method)
  
 # Validate n.cores parameter
  if (!is.numeric(n.cores) || length(n.cores) != 1 || n.cores < 1) {
    stop("n.cores must be a positive integer")
  }
  n.cores <- as.integer(n.cores)
  available_cores <- tryCatch(
    parallel::detectCores(logical = FALSE),
    error = function(e) NA_integer_
  )
  if (!is.na(available_cores) && n.cores > available_cores) {
    warning("n.cores (", n.cores, ") exceeds available physical cores (", 
            available_cores, "). Using ", available_cores, " cores.", call. = FALSE)
    n.cores <- available_cores
  }
  
  if (!inherits(gdata, "GBLUPData")) {
    stop("gdata must be a GBLUPData object")
  }
  
  # Get data
  pheno <- gdata@phenotypes
  M <- gdata@genotypes
  
  if (is.null(M) || nrow(M) == 0) {
    stop("No genotype data available")
  }
  
  # Filter markers by MAF
  M_filtered <- filter_markers_maf(M, min.MAF = min.MAF)
  n_markers <- ncol(M_filtered)
  
  if (verbose) {
    cat("Starting GWAS with", n_markers, "markers\n")
    cat("Method:", method, "\n")
    cat("P3D:", P3D, "\n")
    cat("Using C++:", use.cpp, "\n")
  }
  
  # Align individuals
  common_ids <- intersect(rownames(M_filtered), unique(pheno$GID))
  if (length(common_ids) == 0) {
    stop("No common individuals between genotypes and phenotypes")
  }
  
  M_filtered <- M_filtered[common_ids, , drop = FALSE]
  pheno <- pheno[pheno$GID %in% common_ids, ]
  
  # Aggregate phenotypes by individual (mean across environments)
  pheno_agg <- aggregate(trait ~ GID, data = pheno, FUN = mean, na.rm = TRUE)
  pheno_agg <- pheno_agg[match(common_ids, pheno_agg$GID), ]
  
  # Add principal components if requested
  if (n.PC > 0) {
    if (verbose) cat("Calculating", n.PC, "principal components\n")
    pca <- prcomp(M_filtered, center = TRUE, scale. = TRUE)
    pcs <- pca$x[, seq_len(min(n.PC, ncol(pca$x))), drop = FALSE]
    colnames(pcs) <- paste0("PC", seq_len(ncol(pcs)))
  } else {
    pcs <- NULL
  }
  
  # Build design matrices
  n <- nrow(pheno_agg)
  y <- pheno_agg$trait
  
  # Fixed effects: intercept + PCs
  if (!is.null(pcs)) {
    X <- cbind(1, pcs)
  } else {
    X <- matrix(1, n, 1)
  }
  colnames(X)[1] <- "Intercept"
  
  # Get relationship matrix
  if (!is.null(K_matrices) && length(K_matrices) > 0) {
    if ("A" %in% names(K_matrices)) {
      K <- K_matrices$A
    } else {
      K <- K_matrices[[1]]
    }
    K <- K[common_ids, common_ids]
  } else {
    # No kinship matrix - use identity
    K <- diag(n)
    rownames(K) <- colnames(K) <- common_ids
  }
  
  # Use EMMAX-style efficient algorithm
  if (P3D && use.cpp) {
    marker_results <- gwas_emmax(y, X, M_filtered, K, verbose)
  } else if (P3D) {
    marker_results <- gwas_p3d_r(y, X, M_filtered, K, method, verbose)
  } else {
    marker_results <- gwas_full_r(y, X, M_filtered, K, method, verbose)
  }
  
  # Add marker information
  marker_results$marker <- colnames(M_filtered)
  
  # Calculate genomic control lambda
  min_markers_for_lambda <- 100
  if (n_markers >= min_markers_for_lambda) {
    lambda_val <- calculate_lambda(marker_results$p_value)
  } else {
    warning("Too few markers (", n_markers, ") for reliable lambda estimation. ",
            "Minimum recommended: ", min_markers_for_lambda)
    lambda_val <- NA_real_
  }
  
  if (verbose && !is.na(lambda_val)) {
    cat("Genomic inflation factor (lambda):", round(lambda_val, 3), "\n")
  }
  
  # Create GWASResult object
  gwas_result <- new("GWASResult",
                     results = marker_results,
                     lambda = as.numeric(lambda_val),
                     n_markers = as.numeric(n_markers),
                     min_MAF = as.numeric(min.MAF),
                     method = method)
  
  if (verbose) {
    n_sig <- sum(marker_results$p_value < 0.05 / n_markers, na.rm = TRUE)
    cat("Significant markers (Bonferroni):", n_sig, "\n")
  }
  
  return(gwas_result)
}

#' EMMAX-style Fast GWAS Using C++
#'
#' @keywords internal
gwas_emmax <- function(y, X, M, K, verbose = TRUE) {
  
  n <- length(y)
  m <- ncol(M)
  
  if (verbose) cat("Using EMMAX algorithm (C++)\n")
  
  # Estimate variance components using null model
  emma_result <- tryCatch({
    emma_reml_cpp(y, X, K, tol = 1e-6, max_iter = 100L)
  }, error = function(e) {
    if (verbose) cat("C++ EMMA failed, using R fallback\n")
    tryCatch({
      fit_emma_r(y, X, K, tol = 1e-6, max_iter = 100)
    }, error = function(e2) {
      stop("Both C++ and R EMMA failed. C++ error: ", e$message, 
           "; R error: ", e2$message)
    })
  })
  
  lambda <- emma_result$lambda
  U <- emma_result$U
  d <- emma_result$d
  
  if (verbose) {
    cat("Null model: h2 =", round(emma_result$h2, 3), "\n")
    cat("Testing", m, "markers...\n")
  }
  
  # Use C++ fast marker testing
  results_mat <- tryCatch({
    gwas_emmax_cpp(y, X, M, U, d, lambda)
  }, error = function(e) {
    if (verbose) cat("C++ GWAS failed, using R fallback\n")
    tryCatch({
      gwas_emmax_r(y, X, M, U, d, lambda, verbose)
    }, error = function(e2) {
      stop("Both C++ and R GWAS failed. C++ error: ", e$message,
           "; R error: ", e2$message)
    })
  })
  
  # Convert to data frame
  results <- data.frame(
    beta = results_mat[, 1],
    se = results_mat[, 2],
    t_stat = results_mat[, 3],
    p_value = results_mat[, 4],
    stringsAsFactors = FALSE
  )
  
  return(results)
}

#' EMMAX R Fallback
#' @keywords internal
gwas_emmax_r <- function(y, X, M, U, d, lambda, verbose = TRUE) {
  
  n <- length(y)
  m <- ncol(M)
  p <- ncol(X)
  
  # Transform data
  eta <- as.vector(t(U) %*% y)
  xi <- t(U) %*% X
  
  # Inverse weights
  D_inv <- 1 / (d + lambda)
  
  # Results
  results <- matrix(NA, m, 4)
  colnames(results) <- c("beta", "se", "t_stat", "p_value")
  
  # Degrees of freedom
  df <- n - p - 1
  if (df < 1) df <- 1
  
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = m, style = 3)
  }
  
  for (j in seq_len(m)) {
    marker_trans <- as.vector(t(U) %*% M[, j])
    
    xi_full <- cbind(xi, marker_trans)
    
    xiT_D_xi <- t(xi_full) %*% (D_inv * xi_full)
    xiT_D_eta <- t(xi_full) %*% (D_inv * eta)
    
    beta_full <- tryCatch({
      solve(xiT_D_xi, xiT_D_eta)
    }, error = function(e) {
      solve(xiT_D_xi + diag(1e-8, ncol(xi_full)), xiT_D_eta)
    })
    
    beta_marker <- beta_full[p + 1]
    
    resid <- eta - xi_full %*% beta_full
    rss <- sum(resid^2 * D_inv)
    sigma2_hat <- rss / df
    
    cov_beta <- tryCatch({
      solve(xiT_D_xi)
    }, error = function(e) {
      solve(xiT_D_xi + diag(1e-8, ncol(xi_full)))
    })
    
    se_marker <- sqrt(sigma2_hat * cov_beta[p + 1, p + 1])
    
    t_stat <- beta_marker / se_marker
    p_val <- 2 * pt(abs(t_stat), df, lower.tail = FALSE)
    
    results[j, ] <- c(beta_marker, se_marker, t_stat, p_val)
    
    if (verbose && j %% 100 == 0) {
      setTxtProgressBar(pb, j)
    }
  }
  
  if (verbose) close(pb)
  
  return(results)
}

#' P3D GWAS (R version)
#' @keywords internal
gwas_p3d_r <- function(y, X, M, K, method, verbose) {
  
  n <- length(y)
  m <- ncol(M)
  
  if (verbose) cat("Estimating variance components (P3D)\n")
  
  # Estimate variance components
  vc_result <- estimate_variance_components(
    y = y, X = X, K = K,
    method = "REML", max_iter = 50, tol = 1e-4, verbose = FALSE
  )
  
  sigma2_g <- vc_result$sigma2_g
  sigma2_e <- vc_result$sigma2_e
  
  # Compute V inverse
  V <- sigma2_g * K + sigma2_e * diag(n)
  Vi <- safe_inverse(V)
  if (is.null(Vi)) {
    V <- V + diag(1e-4, n)
    Vi <- solve(V)
  }
  
  # Pre-compute
  XtVi <- t(X) %*% Vi
  XtViX <- XtVi %*% X
  XtViX_inv <- solve(XtViX)
  
  if (verbose) {
    cat("Testing", m, "markers\n")
    pb <- txtProgressBar(min = 0, max = m, style = 3)
  }
  
  results <- test_markers_vectorized(
    y = y, X = X, M = M, K = K, Vi = Vi,
    sigma2_g = sigma2_g, sigma2_e = sigma2_e,
    method = method, P3D = TRUE, n.cores = 1, verbose = verbose
  )
  
  if (verbose) close(pb)
  
  return(results)
}

#' Full GWAS (re-estimate VC per marker)
#' @keywords internal
gwas_full_r <- function(y, X, M, K, method, verbose) {
  
  m <- ncol(M)
  
  if (verbose) {
    cat("Full model (re-estimate VC per marker)\n")
    cat("Testing", m, "markers\n")
  }
  
  results <- test_markers_vectorized(
    y = y, X = X, M = M, K = K, Vi = NULL,
    sigma2_g = NULL, sigma2_e = NULL,
    method = method, P3D = FALSE, n.cores = 1, verbose = verbose
  )
  
  return(results)
}

#' Vectorized Marker Testing
#' @keywords internal
test_markers_vectorized <- function(y, X, M, K, Vi, sigma2_g, sigma2_e,
                                    method, P3D, n.cores, verbose) {
  
  n <- length(y)
  n_markers <- ncol(M)
  p <- ncol(X)
  
  results <- data.frame(
    beta = numeric(n_markers),
    se = numeric(n_markers),
    t_stat = numeric(n_markers),
    p_value = numeric(n_markers),
    stringsAsFactors = FALSE
  )
  
  if (P3D && !is.null(Vi)) {
    XtVi <- t(X) %*% Vi
    XtViX <- XtVi %*% X
    XtViX_inv <- tryCatch(solve(XtViX), error = function(e) {
      solve(XtViX + diag(1e-8, ncol(XtViX)))
    })
  }
  
  if (verbose) {
    pb <- txtProgressBar(min = 0, max = n_markers, style = 3)
  }
  
  for (i in seq_len(n_markers)) {
    marker_geno <- M[, i]
    
    test_result <- tryCatch({
      if (P3D && !is.null(Vi)) {
        test_marker_p3d(y, X, marker_geno, Vi, n)
      } else if (!is.null(K)) {
        test_marker_mixed(y, X, marker_geno, K)
      } else {
        test_marker_ols(y, X, marker_geno)
      }
    }, error = function(e) {
      list(beta = NA, se = NA, t_stat = NA, p_value = NA)
    })
    
    results$beta[i] <- test_result$beta
    results$se[i] <- test_result$se
    results$t_stat[i] <- test_result$t_stat
    results$p_value[i] <- test_result$p_value
    
    if (verbose && i %% 100 == 0) {
      setTxtProgressBar(pb, i)
    }
  }
  
  if (verbose) close(pb)
  
  return(results)
}

#' Test Single Marker with P3D
#' @keywords internal
test_marker_p3d <- function(y, X, marker, Vi, n) {
  
  complete <- complete.cases(y, marker)
  if (sum(complete) < 5) {
    return(list(beta = NA, se = NA, t_stat = NA, p_value = NA))
  }
  
  if (var(marker[complete]) < 1e-10) {
    return(list(beta = NA, se = NA, t_stat = NA, p_value = NA))
  }
  
  y_c <- y[complete]
  X_c <- X[complete, , drop = FALSE]
  marker_c <- marker[complete]
  Vi_c <- Vi[complete, complete]
  
  X_full <- cbind(X_c, marker_c)
  
  XtViX_full <- crossprod(X_full, Vi_c %*% X_full)
  XtViy_full <- crossprod(X_full, Vi_c %*% y_c)
  
  XtViX_inv_full <- tryCatch({
    solve(XtViX_full)
  }, error = function(e) NULL)
  
  if (is.null(XtViX_inv_full)) {
    return(list(beta = NA, se = NA, t_stat = NA, p_value = NA))
  }
  
  beta <- XtViX_inv_full %*% XtViy_full
  se <- sqrt(diag(XtViX_inv_full))
  
  marker_idx <- ncol(X_full)
  marker_beta <- beta[marker_idx]
  marker_se <- se[marker_idx]
  
  if (marker_se == 0 || is.na(marker_se) || !is.finite(marker_se)) {
    return(list(beta = marker_beta, se = marker_se, t_stat = NA, p_value = NA))
  }
  
  t_stat <- marker_beta / marker_se
  df <- sum(complete) - ncol(X_full)
  if (df <= 0) df <- 1
  p_value <- 2 * pt(abs(t_stat), df, lower.tail = FALSE)
  
  return(list(
    beta = as.numeric(marker_beta),
    se = as.numeric(marker_se),
    t_stat = as.numeric(t_stat),
    p_value = as.numeric(p_value)
  ))
}

#' Test Single Marker with OLS
#' @keywords internal
test_marker_ols <- function(y, X, marker) {
  
  complete <- complete.cases(y, marker)
  if (sum(complete) < 5) {
    return(list(beta = NA, se = NA, t_stat = NA, p_value = NA))
  }
  
  if (var(marker[complete]) < 1e-10) {
    return(list(beta = NA, se = NA, t_stat = NA, p_value = NA))
  }
  
  y_c <- y[complete]
  X_c <- cbind(X[complete, , drop = FALSE], marker[complete])
  
  XtX <- crossprod(X_c)
  Xty <- crossprod(X_c, y_c)
  
  XtX_inv <- tryCatch({
    solve(XtX)
  }, error = function(e) NULL)
  
  if (is.null(XtX_inv)) {
    return(list(beta = NA, se = NA, t_stat = NA, p_value = NA))
  }
  
  beta <- XtX_inv %*% Xty
  
  fitted <- X_c %*% beta
  residuals <- y_c - fitted
  df <- sum(complete) - ncol(X_c)
  if (df <= 0) df <- 1
  sigma2 <- sum(residuals^2) / df
  
  se <- sqrt(sigma2 * diag(XtX_inv))
  
  marker_idx <- ncol(X_c)
  marker_beta <- beta[marker_idx]
  marker_se <- se[marker_idx]
  
  if (marker_se == 0 || !is.finite(marker_se)) {
    return(list(beta = marker_beta, se = marker_se, t_stat = NA, p_value = NA))
  }
  
  t_stat <- marker_beta / marker_se
  p_value <- 2 * pt(abs(t_stat), df, lower.tail = FALSE)
  
  return(list(
    beta = as.numeric(marker_beta),
    se = as.numeric(marker_se),
    t_stat = as.numeric(t_stat),
    p_value = as.numeric(p_value)
  ))
}

#' Test Single Marker with Mixed Model
#' @keywords internal
test_marker_mixed <- function(y, X, marker, K) {
  
  complete <- complete.cases(y, marker)
  if (sum(complete) < 5) {
    return(list(beta = NA, se = NA, t_stat = NA, p_value = NA))
  }
  
  y_c <- y[complete]
  X_c <- cbind(X[complete, , drop = FALSE], marker[complete])
  K_c <- K[complete, complete]
  
  vc <- tryCatch({
    estimate_variance_components(y_c, X_c, K_c, max_iter = 20, verbose = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(vc)) {
    return(list(beta = NA, se = NA, t_stat = NA, p_value = NA))
  }
  
  n_effects <- ncol(X_c)
  marker_beta <- vc$beta[n_effects]
  
  V <- vc$sigma2_g * K_c + vc$sigma2_e * diag(nrow(K_c))
  Vi <- tryCatch(solve(V), error = function(e) NULL)
  
  if (is.null(Vi)) {
    return(list(beta = marker_beta, se = NA, t_stat = NA, p_value = NA))
  }
  
  XtViX <- crossprod(X_c, Vi %*% X_c)
  XtViX_inv <- tryCatch(solve(XtViX), error = function(e) NULL)
  
  if (is.null(XtViX_inv)) {
    return(list(beta = marker_beta, se = NA, t_stat = NA, p_value = NA))
  }
  
  marker_se <- sqrt(XtViX_inv[n_effects, n_effects])
  
  if (marker_se == 0 || !is.finite(marker_se)) {
    return(list(beta = marker_beta, se = marker_se, t_stat = NA, p_value = NA))
  }
  
  t_stat <- marker_beta / marker_se
  df <- sum(complete) - n_effects
  if (df <= 0) df <- 1
  p_value <- 2 * pt(abs(t_stat), df, lower.tail = FALSE)
  
  return(list(
    beta = as.numeric(marker_beta),
    se = as.numeric(marker_se),
    t_stat = as.numeric(t_stat),
    p_value = as.numeric(p_value)
  ))
}

#' Estimate Variance Components using REML
#' @keywords internal
estimate_variance_components <- function(y, X, K, method = "REML",
                                         max_iter = 100, tol = 1e-6,
                                         verbose = FALSE) {
  
  n <- length(y)
  p <- ncol(X)
  
  var_y <- var(y, na.rm = TRUE)
  if (var_y == 0) var_y <- 1
  sigma2_g <- var_y * 0.5
  sigma2_e <- var_y * 0.5
  
  eig_K <- eigen(K, symmetric = TRUE)
  U <- eig_K$vectors
  d <- eig_K$values
  d[d < 1e-10] <- 1e-10
  
  Uty <- crossprod(U, y)
  UtX <- crossprod(U, X)
  
  for (iter in seq_len(max_iter)) {
    
    lambda <- sigma2_e / sigma2_g
    D_inv <- 1 / (d + lambda)
    
    XtVinvX <- crossprod(UtX, D_inv * UtX)
    XtVinvX_inv <- tryCatch({
      solve(XtVinvX)
    }, error = function(e) {
      solve(XtVinvX + diag(1e-8, ncol(XtVinvX)))
    })
    
    XtVinvy <- crossprod(UtX, D_inv * Uty)
    beta <- XtVinvX_inv %*% XtVinvy
    
    r <- Uty - UtX %*% beta
    Py <- D_inv * r
    
    yPKPy <- sum(d * Py^2)
    trPK <- sum(d * D_inv) - sum(diag(XtVinvX_inv %*% crossprod(UtX, d * D_inv * UtX)))
    
    yPPy <- sum(Py^2)
    trP <- sum(D_inv) - sum(diag(XtVinvX_inv %*% crossprod(UtX, D_inv * UtX)))
    
    sigma2_g_new <- as.numeric(yPKPy / max(trPK, 1e-10))
    sigma2_e_new <- as.numeric(yPPy / max(trP, 1e-10))
    
    sigma2_g_new <- max(sigma2_g_new, 1e-10)
    sigma2_e_new <- max(sigma2_e_new, 1e-10)
    
    delta_g <- abs(sigma2_g_new - sigma2_g) / max(sigma2_g, 1e-10)
    delta_e <- abs(sigma2_e_new - sigma2_e) / max(sigma2_e, 1e-10)
    
    sigma2_g <- sigma2_g_new
    sigma2_e <- sigma2_e_new
    
    if (max(delta_g, delta_e) < tol) break
  }
  
  h2 <- sigma2_g / (sigma2_g + sigma2_e)
  
  return(list(
    sigma2_g = sigma2_g,
    sigma2_e = sigma2_e,
    h2 = h2,
    beta = as.vector(beta),
    converged = iter < max_iter,
    iterations = iter
  ))
}

#' Calculate Genomic Inflation Factor (Lambda)
#' @keywords internal
calculate_lambda <- function(p_values) {
  
  p_values <- p_values[!is.na(p_values) & p_values > 0 & p_values < 1]
  
  if (length(p_values) < 10) {
    warning("Too few valid p-values (", length(p_values), ") for reliable lambda estimation")
    return(NA)
  }
  
  chi2 <- qchisq(p_values, df = 1, lower.tail = FALSE)
  expected_median <- qchisq(0.5, df = 1)
  lambda <- median(chi2, na.rm = TRUE) / expected_median
  
  return(lambda)
}

#' Manhattan Plot for GWAS Results
#'
#' Creates a publication-quality Manhattan plot with support for single or
#' multiple traits, chromosome coloring, and significance thresholds.
#'
#' @param gwas_result GWASResult object or list of GWASResult objects for multi-trait plot
#' @param chr_info Optional data frame with columns: marker, chr, pos
#' @param threshold Custom -log10(p) threshold line (default: Bonferroni)
#' @param suggestive_threshold Custom suggestive -log10(p) threshold
#' @param significance Significance level for Bonferroni correction (default: 0.05)
#' @param trait_colors Character vector of colors for multi-trait plots
#' @param density_colors Character vector of colors for density gradient (optional)
#' @param point_size Point size (default: 1.2)
#' @param alpha Point transparency (default: 0.7)
#' @param main Plot title
#' @param show_threshold Logical, show threshold lines (default: TRUE)
#' @param ... Additional arguments for plot
#'
#' @return Invisible NULL
#' @export
#' @examples
#' \dontrun{
#' # Single trait
#' plot_manhattan(gwas_result, chr_info = map_data)
#'
#' # Multi-trait
#' gwas_list <- list(Yield = gwas1, Height = gwas2)
#' plot_manhattan(gwas_list, chr_info = map_data)
#' }
plot_manhattan <- function(gwas_result, chr_info = NULL,
                           threshold = NULL, suggestive_threshold = NULL,
                           significance = 0.05,
                           trait_colors = NULL, density_colors = NULL,
                           point_size = 1.2, alpha = 0.7,
                           main = "Manhattan Plot", 
                           show_threshold = TRUE,
                           ...) {
  
 # Handle multiple GWAS results (list) or single result
  if (is.list(gwas_result) && !inherits(gwas_result, "GWASResult")) {
    # Multi-trait Manhattan plot
    return(plot_manhattan_multi_internal(gwas_result, chr_info, threshold,
                                         suggestive_threshold, significance,
                                         trait_colors, point_size, 
                                         alpha, main, show_threshold, ...))
  }
  
  if (!inherits(gwas_result, "GWASResult")) {
    stop("gwas_result must be a GWASResult object or list of GWASResult objects")
  }
  
  results <- gwas_result@results
  results$log10p <- -log10(results$p_value)
  
  # Handle infinite values
  max_finite <- max(results$log10p[is.finite(results$log10p)], na.rm = TRUE)
  results$log10p[is.infinite(results$log10p)] <- max_finite + 1
  results$log10p[is.na(results$log10p)] <- 0
  
  # Calculate threshold
  n_markers <- gwas_result@n_markers
  if (is.null(threshold)) {
    bonferroni <- -log10(significance / n_markers)
  } else {
    bonferroni <- threshold
  }
  
  # Suggestive threshold (default: 1e-5 or 1/n_markers)
  if (is.null(suggestive_threshold)) {
    suggestive_log10 <- -log10(1 / n_markers)
  } else {
    suggestive_log10 <- suggestive_threshold
  }
  
  # Determine y-axis range to include threshold line
  y_max <- max(c(results$log10p, bonferroni * 1.1), na.rm = TRUE)
  
  # Try to extract chromosome info from results if not provided
  has_chr <- FALSE
  if (is.null(chr_info)) {
    # Check if results have chr/CHR column
    chr_col <- intersect(names(results), c("chr", "CHR", "Chr", "chromosome", "Chromosome"))
    pos_col <- intersect(names(results), c("pos", "POS", "Pos", "position", "Position", "bp", "BP"))
    
    if (length(chr_col) > 0 && length(pos_col) > 0) {
      results$chr <- results[[chr_col[1]]]
      results$pos <- results[[pos_col[1]]]
      has_chr <- TRUE
    }
  } else {
    # Merge with provided chr_info
    marker_col <- intersect(names(chr_info), c("marker", "Marker", "SNP", "snp", "ID", "id"))
    if (length(marker_col) > 0) {
      names(chr_info)[names(chr_info) == marker_col[1]] <- "marker"
    }
    results <- merge(results, chr_info, by = "marker", all.x = TRUE)
    has_chr <- "chr" %in% names(results) || "CHR" %in% names(results)
    if (has_chr) {
      chr_col <- intersect(names(results), c("chr", "CHR", "Chr"))
      pos_col <- intersect(names(results), c("pos", "POS", "Pos", "bp", "BP"))
      if (length(chr_col) > 0) results$chr <- results[[chr_col[1]]]
      if (length(pos_col) > 0) results$pos <- results[[pos_col[1]]]
    }
  }
  
  if (has_chr && !all(is.na(results$chr))) {
    # ===== CHROMOSOME-BASED PLOT =====
    results <- results[order(results$chr, results$pos), ]
    
    # Parse chromosome for various naming styles (1A, 1B, 1D, chr1, etc.)
    chr_clean <- gsub("^chr", "", results$chr, ignore.case = TRUE)
    chr_levels <- unique(chr_clean[!is.na(chr_clean)])
    
    # Sort chromosomes naturally
    chr_order <- function(x) {
      num <- as.numeric(gsub("[^0-9]", "", x))
      num[is.na(num)] <- 999
      letter <- gsub("[0-9]", "", x)
      letter_val <- match(toupper(letter), c("A", "B", "C", "D", ""))
      letter_val[is.na(letter_val)] <- 99
      num * 10 + letter_val
    }
    chr_levels <- chr_levels[order(sapply(chr_levels, chr_order))]
    
    results$chr_factor <- factor(chr_clean, levels = chr_levels)
    results <- results[!is.na(results$chr_factor), ]
    
    # Calculate cumulative positions
    chr_info_list <- lapply(split(results, results$chr_factor), function(df) {
      if (nrow(df) == 0) return(NULL)
      data.frame(chr = df$chr_factor[1], chr_len = max(df$pos, na.rm = TRUE))
    })
    chr_info_df <- do.call(rbind, chr_info_list[!sapply(chr_info_list, is.null)])
    chr_info_df <- chr_info_df[match(chr_levels, chr_info_df$chr), ]
    chr_info_df <- chr_info_df[!is.na(chr_info_df$chr), ]
    
    chr_info_df$tot <- c(0, cumsum(as.numeric(chr_info_df$chr_len[-nrow(chr_info_df)])))
    chr_info_df$center <- chr_info_df$tot + chr_info_df$chr_len / 2
    
    results$tot <- chr_info_df$tot[match(results$chr_factor, chr_info_df$chr)]
    results$plot_pos <- results$pos + results$tot
    
    # Alternating chromosome colors
    chr_colors <- rep(c("#4393C3", "#2166AC"), length.out = nrow(chr_info_df))
    results$point_color <- chr_colors[as.numeric(results$chr_factor)]
    
    # Set up plot
    par(mar = c(5, 4, 4, 2) + 0.1)
    
    plot(results$plot_pos, results$log10p,
         xlab = "Chromosome", 
         ylab = expression(-log[10](italic(p))),
         main = main,
         pch = 16, 
         cex = point_size,
         col = adjustcolor(results$point_color, alpha.f = alpha),
         xaxt = "n",
         ylim = c(0, y_max),
         las = 1,
         bty = "l",
         ...)
    
    # Add chromosome labels
    axis(1, at = chr_info_df$center, labels = chr_info_df$chr, tick = FALSE, cex.axis = 0.9)
    
    # Add chromosome color bar at bottom
    usr <- par("usr")
    bar_height <- (usr[4] - usr[3]) * 0.015
    chr_bar_colors <- rep(c("forestgreen", "darkorange"), length.out = nrow(chr_info_df))
    
    for (i in seq_len(nrow(chr_info_df))) {
      x_start <- chr_info_df$tot[i]
      x_end <- chr_info_df$tot[i] + chr_info_df$chr_len[i]
      rect(x_start, usr[3], x_end, usr[3] + bar_height, 
           col = chr_bar_colors[i], border = NA, xpd = FALSE)
    }
    
  } else {
    # ===== INDEX-BASED PLOT (no chromosome info) =====
    # Still make it look good with alternating color blocks
    
    n <- nrow(results)
    block_size <- ceiling(n / 10)  # Create ~10 blocks for visual separation
    results$block <- ceiling(seq_len(n) / block_size)
    block_colors <- rep(c("#4393C3", "#2166AC"), length.out = max(results$block))
    results$point_color <- block_colors[results$block]
    
    par(mar = c(5, 4, 4, 2) + 0.1)
    
    plot(seq_len(n), results$log10p,
         xlab = "Marker Index", 
         ylab = expression(-log[10](italic(p))),
         main = main,
         pch = 16, 
         cex = point_size,
         col = adjustcolor(results$point_color, alpha.f = alpha),
         ylim = c(0, y_max),
         las = 1,
         bty = "l",
         ...)
  }
  
  # Add threshold lines
  if (show_threshold) {
    # Bonferroni/genome-wide significance (red solid line)
    abline(h = bonferroni, col = "#D73027", lty = 1, lwd = 1.5)
    
    # Suggestive threshold (blue dashed line) - only if different from Bonferroni
    if (abs(suggestive_log10 - bonferroni) > 0.5) {
      abline(h = suggestive_log10, col = "#4575B4", lty = 2, lwd = 1)
    }
  }
  
  invisible(NULL)
}

#' Internal function for multi-trait Manhattan plot
#' @keywords internal
plot_manhattan_multi_internal <- function(gwas_list, chr_info = NULL,
                                          threshold = NULL, suggestive_threshold = NULL,
                                          significance = 0.05,
                                          trait_colors = NULL, point_size = 1.2,
                                          alpha = 0.7, main = NULL, 
                                          show_threshold = TRUE, ...) {
  
  trait_names <- names(gwas_list)
  if (is.null(trait_names)) {
    trait_names <- paste0("Trait", seq_along(gwas_list))
  }
  
  # Default colors - vibrant, distinguishable
  if (is.null(trait_colors)) {
    trait_colors <- c("#8B0000", "#0000CD", "#800080", "#DAA520", 
                      "#006400", "#FF4500", "#4B0082", "#2F4F4F")[seq_along(gwas_list)]
  }
  
  # Get n_markers for threshold calculation
  n_markers <- if (inherits(gwas_list[[1]], "GWASResult")) {
    gwas_list[[1]]@n_markers
  } else {
    nrow(gwas_list[[1]])
  }
  
  # Calculate thresholds
  if (is.null(threshold)) {
    bonferroni <- -log10(significance / n_markers)
  } else {
    bonferroni <- threshold
  }
  
  if (is.null(suggestive_threshold)) {
    suggestive_log10 <- -log10(1 / n_markers)
  } else {
    suggestive_log10 <- suggestive_threshold
  }
  
  # Combine all results
  combined <- do.call(rbind, lapply(seq_along(gwas_list), function(i) {
    res <- gwas_list[[i]]
    if (inherits(res, "GWASResult")) {
      df <- res@results
    } else {
      df <- res
    }
    df$trait <- trait_names[i]
    df$log10p <- -log10(df$p_value)
    df$log10p[!is.finite(df$log10p)] <- NA
    df
  }))
  
  combined <- combined[!is.na(combined$log10p), ]
  
  # Y-axis range
  y_max <- max(c(combined$log10p, bonferroni * 1.1), na.rm = TRUE)
  
  # Check for chromosome info in data
  has_chr <- FALSE
  if (!is.null(chr_info)) {
    marker_col <- intersect(names(chr_info), c("marker", "Marker", "SNP", "snp"))
    if (length(marker_col) > 0) names(chr_info)[names(chr_info) == marker_col[1]] <- "marker"
    combined <- merge(combined, chr_info, by = "marker", all.x = TRUE)
    has_chr <- TRUE
  } else {
    chr_col <- intersect(names(combined), c("chr", "CHR", "Chr"))
    pos_col <- intersect(names(combined), c("pos", "POS", "Pos", "bp", "BP"))
    if (length(chr_col) > 0 && length(pos_col) > 0) {
      combined$chr <- combined[[chr_col[1]]]
      combined$pos <- combined[[pos_col[1]]]
      has_chr <- TRUE
    }
  }
  
  if (has_chr && "chr" %in% names(combined) && !all(is.na(combined$chr))) {
    # ===== CHROMOSOME-BASED MULTI-TRAIT PLOT =====
    combined <- combined[order(combined$chr, combined$pos), ]
    
    chr_clean <- gsub("^chr", "", combined$chr, ignore.case = TRUE)
    chr_levels <- unique(chr_clean[!is.na(chr_clean)])
    
    chr_order <- function(x) {
      num <- as.numeric(gsub("[^0-9]", "", x))
      num[is.na(num)] <- 999
      letter <- gsub("[0-9]", "", x)
      letter_val <- match(toupper(letter), c("A", "B", "C", "D", ""))
      letter_val[is.na(letter_val)] <- 99
      num * 10 + letter_val
    }
    chr_levels <- chr_levels[order(sapply(chr_levels, chr_order))]
    
    combined$chr_factor <- factor(chr_clean, levels = chr_levels)
    combined <- combined[!is.na(combined$chr_factor), ]
    
    chr_info_list <- lapply(split(combined, combined$chr_factor), function(df) {
      if (nrow(df) == 0) return(NULL)
      data.frame(chr = df$chr_factor[1], chr_len = max(df$pos, na.rm = TRUE))
    })
    chr_info_df <- do.call(rbind, chr_info_list[!sapply(chr_info_list, is.null)])
    chr_info_df <- chr_info_df[match(chr_levels, chr_info_df$chr), ]
    chr_info_df <- chr_info_df[!is.na(chr_info_df$chr), ]
    
    chr_info_df$tot <- c(0, cumsum(as.numeric(chr_info_df$chr_len[-nrow(chr_info_df)])))
    chr_info_df$center <- chr_info_df$tot + chr_info_df$chr_len / 2
    
    combined$tot <- chr_info_df$tot[match(combined$chr_factor, chr_info_df$chr)]
    combined$plot_pos <- combined$pos + combined$tot
    
    x_var <- combined$plot_pos
    x_range <- range(x_var, na.rm = TRUE)
    x_labels <- chr_info_df$chr
    x_at <- chr_info_df$center
    
  } else {
    # ===== INDEX-BASED MULTI-TRAIT PLOT =====
    combined$plot_pos <- seq_len(nrow(combined))
    x_var <- combined$plot_pos
    x_range <- range(x_var, na.rm = TRUE)
    x_labels <- NULL
    x_at <- NULL
    chr_info_df <- NULL
  }
  
  # Set up plot
  if (is.null(main)) main <- "Multi-Trait GWAS"
  
  par(mar = c(5, 4, 4, 2) + 0.1)
  
  plot(x_range, c(0, y_max), type = "n",
       xlab = if (!is.null(x_labels)) "Chromosome" else "Marker Index",
       ylab = expression(-log[10](italic(p))),
       main = main,
       xaxt = "n",
       las = 1,
       bty = "l",
       ...)
  
  # Plot points for each trait
  for (i in seq_along(trait_names)) {
    trait_data <- combined[combined$trait == trait_names[i], ]
    points(trait_data$plot_pos, trait_data$log10p,
           pch = 16, cex = point_size,
           col = adjustcolor(trait_colors[i], alpha.f = alpha))
  }
  
  # Add chromosome axis and color bar
  if (!is.null(x_labels)) {
    axis(1, at = x_at, labels = x_labels, tick = FALSE, cex.axis = 0.9)
    
    # Chromosome color bar
    usr <- par("usr")
    bar_height <- (usr[4] - usr[3]) * 0.015
    chr_bar_colors <- rep(c("forestgreen", "darkorange"), length.out = nrow(chr_info_df))
    
    for (j in seq_len(nrow(chr_info_df))) {
      x_start <- chr_info_df$tot[j]
      x_end <- chr_info_df$tot[j] + chr_info_df$chr_len[j]
      rect(x_start, usr[3], x_end, usr[3] + bar_height,
           col = chr_bar_colors[j], border = NA, xpd = FALSE)
    }
  } else {
    axis(1)
  }
  
  # Threshold lines
  if (show_threshold) {
    abline(h = bonferroni, col = "#D73027", lty = 1, lwd = 1.5)
    if (abs(suggestive_log10 - bonferroni) > 0.5) {
      abline(h = suggestive_log10, col = "#4575B4", lty = 2, lwd = 1)
    }
  }
  
  # Legend at top
  legend("top", legend = trait_names, col = trait_colors[seq_along(trait_names)],
         pch = 16, bty = "n", horiz = TRUE, cex = 0.9, pt.cex = 1.2)
  
  invisible(NULL)
}

#' QQ Plot for GWAS Results
#'
#' @param gwas_result GWASResult object
#' @param main Plot title
#' @param ... Additional arguments for plot
#'
#' @return Invisible NULL
#' @export
plot_qq <- function(gwas_result, main = "QQ Plot", ...) {
  
  if (!inherits(gwas_result, "GWASResult")) {
    stop("gwas_result must be a GWASResult object")
  }
  
  p_obs <- gwas_result@results$p_value
  p_obs <- p_obs[!is.na(p_obs) & p_obs > 0 & p_obs < 1]
  
  if (length(p_obs) == 0) {
    stop("No valid p-values to plot")
  }
  
  p_obs <- sort(p_obs)
  n <- length(p_obs)
  p_exp <- (seq_len(n)) / (n + 1)
  
  log10_obs <- -log10(p_obs)
  log10_exp <- -log10(p_exp)
  
  plot(log10_exp, log10_obs,
       xlab = "Expected -log10(p-value)",
       ylab = "Observed -log10(p-value)",
       main = main,
       pch = 16, col = rgb(0, 0, 1, 0.5), ...)
  
  abline(0, 1, col = "red", lwd = 2)
  
  lambda <- gwas_result@lambda
  if (!is.na(lambda)) {
    legend("topleft", 
           legend = paste("lambda =", round(lambda, 3)),
           bty = "n")
  }
  
  add_qq_confidence_band(n)
  
  invisible(NULL)
}

#' Add Confidence Band to QQ Plot
#' @keywords internal
add_qq_confidence_band <- function(n, conf_level = 0.95) {
  
  alpha <- 1 - conf_level
  p_exp <- (seq_len(n)) / (n + 1)
  
  lower <- qbeta(alpha / 2, seq_len(n), n:1 + 1)
  upper <- qbeta(1 - alpha / 2, seq_len(n), n:1 + 1)
  
  log10_exp <- -log10(p_exp)
  log10_lower <- -log10(upper)
  log10_upper <- -log10(lower)
  
  polygon(c(log10_exp, rev(log10_exp)),
          c(log10_lower, rev(log10_upper)),
          col = rgb(0.5, 0.5, 0.5, 0.2),
          border = NA)
}

#' Get Top Markers from GWAS
#'
#' @param gwas_result GWASResult object
#' @param n Number of top markers to return
#' @param threshold Optional p-value threshold
#'
#' @return Data frame with top markers
#' @export
get_top_markers <- function(gwas_result, n = 20, threshold = NULL) {
  
  if (!inherits(gwas_result, "GWASResult")) {
    stop("gwas_result must be a GWASResult object")
  }
  
  results <- gwas_result@results
  results <- results[order(results$p_value), ]
  
  if (!is.null(threshold)) {
    results <- results[results$p_value < threshold, ]
  }
  
  if (nrow(results) > n) {
    results <- results[seq_len(n), ]
  }
  
  results$log10p <- -log10(results$p_value)
  bonf_threshold <- 0.05 / gwas_result@n_markers
  results$significant <- results$p_value < bonf_threshold
  
  return(results)
}

#' Get Significant Markers
#'
#' @param gwas_result GWASResult object
#' @param method Correction method: "bonferroni", "fdr", or "none"
#' @param alpha Significance level
#'
#' @return Data frame with significant markers
#' @export
get_significant_markers <- function(gwas_result, 
                                    method = c("bonferroni", "fdr", "none"), 
                                    alpha = 0.05) {
  
  method <- match.arg(method)
  
  if (!inherits(gwas_result, "GWASResult")) {
    stop("gwas_result must be a GWASResult object")
  }
  
  results <- gwas_result@results
  
  if (method == "bonferroni") {
    threshold <- alpha / gwas_result@n_markers
    sig_markers <- results[results$p_value < threshold, ]
    
  } else if (method == "fdr") {
    # Correct Benjamini-Hochberg FDR procedure
    # Find the largest k such that P(k) <= k/m * alpha
    p_sorted <- sort(results$p_value)
    m <- length(p_sorted)
    k <- seq_len(m)
    
    # Element-wise comparison: P(k) <= k/m * alpha
    significant <- p_sorted <= (k / m) * alpha
    
    # Find largest k where this holds
    if (any(significant)) {
      max_k <- max(which(significant))
      fdr_threshold <- p_sorted[max_k]
      sig_markers <- results[results$p_value <= fdr_threshold, ]
    } else {
      # No significant markers
      sig_markers <- results[0, , drop = FALSE]
    }
    
  } else if (method == "none") {
    sig_markers <- results[results$p_value < alpha, ]
  }
  
  if (nrow(sig_markers) > 0) {
    sig_markers <- sig_markers[order(sig_markers$p_value), ]
    sig_markers$log10p <- -log10(sig_markers$p_value)
  }
  
  return(sig_markers)
}

#' Apply Genomic Control Correction
#'
#' @param gwas_result GWASResult object
#'
#' @return GWASResult object with corrected p-values
#' @export
apply_genomic_control <- function(gwas_result) {
  
  if (!inherits(gwas_result, "GWASResult")) {
    stop("gwas_result must be a GWASResult object")
  }
  
  lambda <- gwas_result@lambda
  
  if (is.na(lambda) || lambda <= 0) {
    warning("Invalid lambda, cannot apply genomic control")
    return(gwas_result)
  }
  
  chi2 <- qchisq(gwas_result@results$p_value, df = 1, lower.tail = FALSE)
  chi2_corrected <- chi2 / lambda
  p_corrected <- pchisq(chi2_corrected, df = 1, lower.tail = FALSE)
  
  new_result <- gwas_result
  new_result@results$p_value_original <- gwas_result@results$p_value
  new_result@results$p_value <- p_corrected
  new_result@lambda <- 1.0
  
  return(new_result)
}
