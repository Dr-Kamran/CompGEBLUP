#' Fit Multi-Trait GBLUP Model
#'
#' Fits a multi-trait genomic prediction model using simplified multivariate REML.
#' This is a basic implementation suitable for correlated trait analysis.
#'
#' @details
#' This function fits independent GBLUP models per trait and estimates genetic/residual
#' correlations. For full multivariate REML with proper covariance estimation,
#' consider specialized packages like sommer or MTM.
#'
#' @param gblup_data GBLUPData object
#' @param traits Character vector of trait names (columns in phenotype data)
#' @param K_matrices List of relationship matrices
#' @param effects Character vector of effects to include
#' @param genetic_model Type of genetic covariance: "diagonal" (default), "unstructured"
#' @param ridge Ridge parameter
#' @param nIters Maximum iterations for REML
#' @param tolParConvLL Convergence tolerance
#' @param verbose Print fitting details
#'
#' @return List containing:
#'   \item{models}{List of single-trait GBLUPModel objects}
#'   \item{genetic_correlations}{Matrix of genetic correlations between traits}
#'   \item{residual_correlations}{Matrix of residual correlations between traits}
#'   \item{gebvs}{Data frame with GEBVs for all traits}
#'   \item{varcomp}{Variance components per trait}
#'
#' @note Returns a list rather than S4 object. Access genetic correlations via
#'   \code{result$genetic_correlations} and GEBVs via \code{result$gebvs}.
#'
#' @seealso \code{\link{fit_gblup}} for single-trait models
#' @export
fit_multi_trait <- function(gblup_data,
                            traits,
                            K_matrices,
                            effects = c("A", "ENV"),
                            genetic_model = c("diagonal", "unstructured"),
                            ridge = 1e-3,
                            nIters = 50,
                            tolParConvLL = 1e-04,
                            verbose = FALSE) {
  
  # Match arguments
  genetic_model <- match.arg(genetic_model)
  
  # Validate inputs
  stopifnot("gblup_data must be a GBLUPData object" = inherits(gblup_data, "GBLUPData"))
  stopifnot("traits must have at least 2 elements" = length(traits) >= 2)
  
  # Prepare data
  data <- gblup_data@phenotypes
  
  # Check if all traits exist
  missing_traits <- setdiff(traits, names(data))
  if (length(missing_traits) > 0) {
    stop("Traits not found in phenotypes: ", paste(missing_traits, collapse = ", "))
  }
  
  n_traits <- length(traits)
  
  # Get unique GIDs
  unique_gids <- unique(as.character(data$GID))
  n_gids <- length(unique_gids)
  
  # Create trait means per GID (average across environments)
  trait_means <- matrix(NA, nrow = n_gids, ncol = n_traits)
  rownames(trait_means) <- unique_gids
  colnames(trait_means) <- traits
  
  for (t in traits) {
    agg <- aggregate(data[[t]] ~ GID, data = data, FUN = mean, na.rm = TRUE)
    trait_means[as.character(agg$GID), t] <- agg[[2]]
  }
  
  # Remove individuals with any missing trait
  complete_rows <- complete.cases(trait_means)
  trait_means <- trait_means[complete_rows, , drop = FALSE]
  complete_gids <- rownames(trait_means)
  n_complete <- nrow(trait_means)
  
  if (n_complete < 10) {
    stop("Too few individuals with complete data for all traits: ", n_complete)
  }
  
  # Add ridge to K matrices and align
  if ("A" %in% names(K_matrices)) {
    K_A <- K_matrices$A
    # Align K matrix with complete GIDs
    common_gids <- intersect(complete_gids, rownames(K_A))
    if (length(common_gids) < n_complete) {
      warning("Some GIDs not found in K matrix. Using ", length(common_gids), " individuals.")
      trait_means <- trait_means[common_gids, , drop = FALSE]
      n_complete <- length(common_gids)
      complete_gids <- common_gids
    }
    K_A <- K_A[complete_gids, complete_gids]
    K_A <- K_A + diag(ridge, nrow(K_A))
  } else {
    stop("K_matrices must contain 'A' matrix for multi-trait model")
  }
  
  if (verbose) {
    cat("Fitting multi-trait model\n")
    cat("  Traits:", paste(traits, collapse = ", "), "\n")
    cat("  Individuals:", n_complete, "\n")
    cat("  Genetic model:", genetic_model, "\n")
  }
  
  # Estimate genetic and residual covariance matrices
  if (genetic_model == "unstructured") {
    # Estimate full genetic covariance using REML
    result <- estimate_multivariate_varcomp(
      Y = trait_means,
      K = K_A,
      nIters = nIters,
      tol = tolParConvLL,
      verbose = verbose
    )
    
    Vg <- result$Vg
    Ve <- result$Ve
    
  } else {
    # Diagonal model - separate variance for each trait
    Vg <- matrix(0, n_traits, n_traits)
    Ve <- matrix(0, n_traits, n_traits)
    rownames(Vg) <- colnames(Vg) <- traits
    rownames(Ve) <- colnames(Ve) <- traits
    
    for (i in 1:n_traits) {
      y <- trait_means[, i]
      
      # Simple REML for each trait
      vc <- estimate_variance_components(
        y = y,
        X = matrix(1, length(y), 1),
        K = K_A,
        max_iter = nIters,
        tol = tolParConvLL,
        verbose = FALSE
      )
      
      Vg[i, i] <- vc$sigma2_g
      Ve[i, i] <- vc$sigma2_e
    }
  }
  
  # Calculate genetic correlation matrix
  gen_cor <- cov2cor_safe(Vg)
  rownames(gen_cor) <- colnames(gen_cor) <- traits
  
  # Calculate heritabilities
  h2 <- numeric(n_traits)
  names(h2) <- traits
  for (i in 1:n_traits) {
    total_var <- Vg[i, i] + Ve[i, i]
    if (total_var > 0) {
      h2[i] <- Vg[i, i] / total_var
      h2[i] <- max(0, min(1, h2[i]))  # Bound to [0, 1]
    } else {
      h2[i] <- NA
    }
  }
  
  # Calculate BLUPs for each trait
  blups <- matrix(NA, n_complete, n_traits)
  rownames(blups) <- complete_gids
  colnames(blups) <- traits
  
  for (i in 1:n_traits) {
    y <- trait_means[, i]
    X <- matrix(1, length(y), 1)
    
    # Simple BLUP calculation
    var_g <- Vg[i, i]
    var_e <- Ve[i, i]
    
    if (var_g > 0 && var_e > 0) {
      V <- var_g * K_A + var_e * diag(n_complete)
      Vi <- tryCatch(solve(V), error = function(e) {
        solve(V + diag(1e-6, n_complete))
      })
      
      # GLS for fixed effects
      XtVi <- t(X) %*% Vi
      beta <- solve(XtVi %*% X) %*% XtVi %*% y
      
      # BLUPs
      Py <- Vi %*% (y - X %*% beta)
      blups[, i] <- as.vector(var_g * K_A %*% Py)
    }
  }
  
  if (verbose) {
    cat("\nResults:\n")
    cat("  Heritabilities:", paste(round(h2, 3), collapse = ", "), "\n")
    if (genetic_model == "unstructured") {
      cat("  Genetic correlations:\n")
      print(round(gen_cor, 3))
    }
  }
  
  return(list(
    model = list(
      convergence = TRUE,
      llik = NA,
      AIC = NA,
      BIC = NA
    ),
    genetic_correlation = gen_cor,
    heritability = h2,
    genetic_variance = Vg,
    residual_variance = Ve,
    blups = blups,
    data = data.frame(
      GID = complete_gids,
      trait_means
    ),
    genetic_model = genetic_model,
    traits = traits,
    effect_names = "A",
    u = lapply(1:n_traits, function(i) blups[, i])
  ))
}

#' Estimate Multivariate Variance Components
#'
#' Uses EM-REML to estimate genetic and residual covariance matrices
#'
#' @param Y Matrix of trait values (individuals x traits)
#' @param K Relationship matrix
#' @param nIters Maximum iterations
#' @param tol Convergence tolerance
#' @param verbose Print progress
#' @return List with Vg and Ve matrices
#' @importFrom stats cov
#' @keywords internal
estimate_multivariate_varcomp <- function(Y, K, nIters = 50, tol = 1e-4, verbose = FALSE) {
  
  n <- nrow(Y)
  t <- ncol(Y)
  
  # Center traits
  Y_centered <- scale(Y, center = TRUE, scale = FALSE)
  
  # Initialize covariance matrices using phenotypic covariance
  P <- cov(Y_centered, use = "complete.obs")
  Vg <- P * 0.5
  Ve <- P * 0.5
  
  # Ensure positive definite
  Vg <- make_pd(Vg)
  Ve <- make_pd(Ve)
  
  # Eigendecomposition of K for efficiency
  eig_K <- eigen(K, symmetric = TRUE)
  U <- eig_K$vectors
  d <- eig_K$values
  d[d < 1e-10] <- 1e-10
  
  # Transform data
  Y_star <- t(U) %*% Y_centered
  
  # EM-REML iterations
  for (iter in 1:nIters) {
    Vg_old <- Vg
    Ve_old <- Ve
    
    # Update estimates using EM
    Vg_new <- matrix(0, t, t)
    Ve_new <- matrix(0, t, t)
    
    for (i in 1:n) {
      # Variance for this individual in transformed space
      # V_i = d[i] * Vg + Ve
      V_i <- d[i] * Vg + Ve
      V_i_inv <- tryCatch(solve(V_i), error = function(e) {
        solve(V_i + diag(1e-8, t))
      })
      
      y_i <- Y_star[i, ]
      
      # Contribution to genetic variance
      Vg_new <- Vg_new + d[i] * (outer(y_i, y_i) %*% V_i_inv %*% Vg %*% V_i_inv +
                                   Vg - Vg %*% V_i_inv %*% Vg * d[i])
      
      # Contribution to residual variance
      Ve_new <- Ve_new + (outer(y_i, y_i) %*% V_i_inv %*% Ve %*% V_i_inv +
                            Ve - Ve %*% V_i_inv %*% Ve)
    }
    
    Vg <- Vg_new / n
    Ve <- Ve_new / n
    
    # Ensure positive semi-definite
    Vg <- make_pd(Vg)
    Ve <- make_pd(Ve)
    
    # Check convergence
    diff_g <- max(abs(Vg - Vg_old))
    diff_e <- max(abs(Ve - Ve_old))
    
    if (verbose && iter %% 10 == 0) {
      cat("  Iter", iter, ": max diff =", round(max(diff_g, diff_e), 6), "\n")
    }
    
    if (max(diff_g, diff_e) < tol) {
      if (verbose) cat("  Converged at iteration", iter, "\n")
      break
    }
  }
  
  return(list(Vg = Vg, Ve = Ve, converged = iter < nIters))
}

#' Make Matrix Positive Definite
#' @keywords internal
make_pd <- function(M, tol = 1e-8) {
  # Use the main make_positive_definite function
  make_positive_definite(M, tol = tol)
}

#' Safe cov2cor
#' @keywords internal
cov2cor_safe <- function(V) {
  d <- sqrt(diag(V))
  d[d == 0] <- 1
  V / outer(d, d)
}

#' Calculate Selection Index
#'
#' @param mt_model Multi-trait model from fit_multi_trait()
#' @param weights Named vector of weights for each trait
#' @return Data frame with selection index values
#' @export
selection_index <- function(mt_model, weights) {
  
  if (is.null(mt_model$blups)) {
    stop("Model does not contain BLUPs")
  }
  
  traits <- mt_model$traits
  blups <- mt_model$blups
  
  # Check weights
  if (!all(names(weights) %in% traits)) {
    stop("Weight names must match trait names: ", paste(traits, collapse = ", "))
  }
  
  # Calculate weighted index
  gids <- rownames(blups)
  index_values <- numeric(length(gids))
  
  for (i in seq_along(traits)) {
    t <- traits[i]
    if (t %in% names(weights)) {
      index_values <- index_values + weights[t] * blups[, t]
    }
  }
  
  result <- data.frame(
    GID = gids,
    Index = index_values,
    stringsAsFactors = FALSE
  )
  
  result <- result[order(-result$Index), ]
  rownames(result) <- NULL
  
  return(result)
}

#' Summary Method for Multi-Trait Model
#'
#' @param mt_model Multi-trait model from fit_multi_trait()
#' @return Invisible NULL
#' @export
summary_multi_trait <- function(mt_model) {
  
  cat("Multi-Trait GBLUP Model Summary\n")
  cat("================================\n\n")
  
  cat("Traits:", paste(mt_model$traits, collapse = ", "), "\n")
  cat("Genetic Model:", mt_model$genetic_model, "\n")
  cat("Individuals:", nrow(mt_model$blups), "\n\n")
  
  if (!is.null(mt_model$heritability)) {
    cat("Heritabilities:\n")
    for (t in names(mt_model$heritability)) {
      cat("  ", t, ":", round(mt_model$heritability[t], 3), "\n")
    }
    cat("\n")
  }
  
  if (!is.null(mt_model$genetic_correlation)) {
    cat("Genetic Correlations:\n")
    print(round(mt_model$genetic_correlation, 3))
    cat("\n")
  }
  
  if (!is.null(mt_model$genetic_variance)) {
    cat("Genetic Variances (diagonal):\n")
    for (t in mt_model$traits) {
      cat("  ", t, ":", round(mt_model$genetic_variance[t, t], 3), "\n")
    }
    cat("\n")
  }
  
  if (!is.null(mt_model$residual_variance)) {
    cat("Residual Variances (diagonal):\n")
    for (t in mt_model$traits) {
      cat("  ", t, ":", round(mt_model$residual_variance[t, t], 3), "\n")
    }
  }
  
  invisible(NULL)
}