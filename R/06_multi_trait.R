#' Fit Multi-Trait GBLUP Model
#'
#' Fits a multi-trait genomic prediction model using multivariate EM-REML
#' with spectral decomposition for efficient computation.
#'
#' @details
#' This function estimates genetic and residual covariance matrices across
#' multiple traits using EM-REML in the spectral space of the kinship matrix.
#' For the unstructured model, it uses proper multivariate EM-REML with
#' step-size damping for stable convergence. Traits on very different scales
#' are automatically standardized (with back-transformation of results).
#'
#' @param gblup_data GBLUPData object
#' @param traits Character vector of trait names (columns in phenotype data)
#' @param K_matrices List of relationship matrices (must contain 'A')
#' @param effects Character vector of effects to include
#' @param genetic_model Type of genetic covariance: "diagonal" (default), "unstructured"
#' @param ridge Ridge parameter for K matrices
#' @param nIters Maximum iterations for EM-REML
#' @param tolParConvLL Convergence tolerance (on relative log-likelihood change)
#' @param verbose Print fitting details
#'
#' @return List containing:
#'   \item{model}{List with convergence info and log-likelihood}
#'   \item{genetic_correlation}{Matrix of genetic correlations between traits}
#'   \item{heritability}{Named numeric vector of per-trait heritabilities}
#'   \item{genetic_variance}{Genetic covariance matrix (Vg)}
#'   \item{residual_variance}{Residual covariance matrix (Ve)}
#'   \item{blups}{Matrix of BLUPs (individuals x traits)}
#'   \item{data}{Data frame with GID and trait means used}
#'   \item{genetic_model}{Model type used}
#'   \item{traits}{Trait names}
#'   \item{effect_names}{Effects included}
#'   \item{u}{List of BLUP vectors per trait}
#'
#' @note Returns a list rather than S4 object. Access genetic correlations via
#'   \code{result$genetic_correlation} and GEBVs via \code{result$blups}.
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

  for (tr in traits) {
    agg <- aggregate(data[[tr]] ~ GID, data = data, FUN = mean, na.rm = TRUE)
    trait_means[as.character(agg$GID), tr] <- agg[[2]]
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

  # --- Auto-standardization for scale-heterogeneous traits ---
  trait_sds <- apply(trait_means, 2, sd, na.rm = TRUE)
  scale_ratio <- max(trait_sds) / max(min(trait_sds), 1e-10)
  did_standardize <- FALSE
  Y_means_orig <- colMeans(trait_means, na.rm = TRUE)
  Y_sds_orig   <- trait_sds

  if (scale_ratio > 10 && genetic_model == "unstructured") {
    if (verbose) {
      cat("  Auto-standardizing traits (scale ratio =", round(scale_ratio, 1), ")\n")
    }
    trait_means <- scale(trait_means)
    did_standardize <- TRUE
  }

  if (verbose) {
    cat("Fitting multi-trait model\n")
    cat("  Traits:", paste(traits, collapse = ", "), "\n")
    cat("  Individuals:", n_complete, "\n")
    cat("  Genetic model:", genetic_model, "\n")
  }

  # Estimate genetic and residual covariance matrices
  if (genetic_model == "unstructured") {
    # Estimate full genetic covariance using EM-REML
    result <- estimate_multivariate_varcomp(
      Y = trait_means,
      K = K_A,
      nIters = nIters,
      tol = tolParConvLL,
      verbose = verbose
    )

    Vg <- result$Vg
    Ve <- result$Ve

    if (!result$converged) {
      warning("EM-REML did not converge after ", nIters, " iterations. ",
              "Results may be unreliable. Consider increasing nIters or ",
              "checking trait scaling.", call. = FALSE)
    }

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

  # --- Back-transform if standardized ---
  if (did_standardize) {
    D_sd <- diag(Y_sds_orig)
    Vg <- D_sd %*% Vg %*% D_sd
    Ve <- D_sd %*% Ve %*% D_sd
    # Restore original data for BLUPs
    for (j in 1:n_traits) {
      trait_means[, j] <- trait_means[, j] * Y_sds_orig[j] + Y_means_orig[j]
    }
  }

  rownames(Vg) <- colnames(Vg) <- traits
  rownames(Ve) <- colnames(Ve) <- traits

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

  # Warn about suspiciously high heritabilities
  if (any(h2 > 0.99, na.rm = TRUE)) {
    high_h2 <- traits[which(h2 > 0.99)]
    warning("Suspiciously high heritability (>0.99) for: ",
            paste(high_h2, collapse = ", "),
            ". This may indicate model mis-specification or EM-REML convergence issues.",
            call. = FALSE)
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

  # Build convergence info
  conv_info <- list(
    convergence = if (genetic_model == "unstructured") result$converged else TRUE,
    llik = if (genetic_model == "unstructured") result$loglik else NA,
    AIC = NA,
    BIC = NA
  )

  return(list(
    model = conv_info,
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

#' Estimate Multivariate Variance Components via EM-REML
#'
#' Uses EM-REML with spectral decomposition and step-size damping to estimate
#' genetic and residual covariance matrices for multi-trait models.
#'
#' @details
#' The algorithm works in the spectral space of the kinship matrix K = UDU'.
#' Transformed data Y* = U'Y yields independent observations with
#' Var(y*_i) = d_i * Vg + Ve.
#'
#' The EM-REML updates follow Meyer (1997) and Lee & van der Werf (2006):
#' \deqn{Vg_{new} = (1/n) \sum_i [Vg - d_i Vg V_i^{-1} Vg + d_i Vg V_i^{-1} y_i y_i' V_i^{-1} Vg]}
#' \deqn{Ve_{new} = (1/n) \sum_i [Ve - Ve V_i^{-1} Ve + Ve V_i^{-1} y_i y_i' V_i^{-1} Ve]}
#'
#' Step-size damping prevents oscillation:
#' \deqn{Vg = \alpha * Vg_{new} + (1-\alpha) * Vg_{old}}
#'
#' @param Y Matrix of trait values (individuals x traits), should be centered/scaled
#' @param K Relationship matrix (n x n)
#' @param nIters Maximum EM iterations (default 500)
#' @param tol Convergence tolerance on relative log-likelihood change (default 1e-6)
#' @param alpha EM damping factor in (0, 1]. Lower = more damping = slower but stable.
#'   Default 0.7 works well for most cases.
#' @param verbose Print progress every 10 iterations
#' @return List with Vg, Ve, converged flag, iteration count, and REML log-likelihood
#' @importFrom stats cov var
#' @keywords internal
estimate_multivariate_varcomp <- function(Y, K, nIters = 500, tol = 1e-6,
                                          alpha = 0.7, verbose = FALSE) {

  n <- nrow(Y)
  n_traits <- ncol(Y)

  # Validate alpha
  alpha <- max(0.1, min(1.0, alpha))

  # Center traits (remove intercept)
  Y_centered <- scale(Y, center = TRUE, scale = FALSE)

  # ---- Eigendecomposition of K for spectral transformation ----
  eig_K <- eigen(K, symmetric = TRUE)
  U <- eig_K$vectors
  d <- eig_K$values
  d[d < 1e-10] <- 1e-10  # floor small eigenvalues

  # Transform data to spectral space: Y* = U' Y
  Y_star <- crossprod(U, Y_centered)  # n x n_traits

  # ---- Initialize Vg and Ve using per-trait grid search ----
  # Much better starting values than naive 50/50 split of P
  Vg <- matrix(0, n_traits, n_traits)
  Ve <- matrix(0, n_traits, n_traits)

  for (j in 1:n_traits) {
    y_j_star <- Y_star[, j]
    vy <- var(Y_centered[, j], na.rm = TRUE)
    if (vy < 1e-10) vy <- 1

    # Quick 1D REML via grid search on h2
    best_ll <- -Inf
    best_h2 <- 0.5
    for (h2_try in seq(0.1, 0.9, by = 0.1)) {
      vg_try <- vy * h2_try
      ve_try <- vy * (1 - h2_try)
      ll <- -0.5 * sum(log(d * vg_try + ve_try) +
                          y_j_star^2 / (d * vg_try + ve_try))
      if (is.finite(ll) && ll > best_ll) {
        best_ll <- ll
        best_h2 <- h2_try
      }
    }
    Vg[j, j] <- vy * best_h2
    Ve[j, j] <- vy * (1 - best_h2)
  }

  # Initialize off-diagonals from phenotypic correlations (conservative)
  P <- cov(Y_centered, use = "complete.obs")
  for (i in 1:(n_traits - 1)) {
    for (j in (i + 1):n_traits) {
      rp <- P[i, j] / sqrt(max(P[i, i], 1e-10) * max(P[j, j], 1e-10))
      rp <- max(-0.95, min(0.95, rp))  # bound to avoid singularity
      Vg[i, j] <- Vg[j, i] <- rp * sqrt(Vg[i, i] * Vg[j, j]) * 0.5
      Ve[i, j] <- Ve[j, i] <- rp * sqrt(Ve[i, i] * Ve[j, j]) * 0.5
    }
  }

  # Ensure positive definite
  Vg <- make_pd(Vg)
  Ve <- make_pd(Ve)

  if (verbose) {
    cat("  EM-REML initialization:\n")
    cat("    Per-trait h2:", paste(round(diag(Vg) / (diag(Vg) + diag(Ve)), 3),
                                   collapse = ", "), "\n")
  }

  # ---- EM-REML iterations ----
  loglik <- -Inf
  converged <- FALSE
  final_iter <- nIters

  for (iter in 1:nIters) {
    Vg_old <- Vg
    Ve_old <- Ve
    loglik_old <- loglik

    # Accumulators for EM updates
    Vg_acc <- matrix(0, n_traits, n_traits)
    Ve_acc <- matrix(0, n_traits, n_traits)
    loglik <- 0

    for (i in 1:n) {
      # Per-individual variance in spectral space: V_i = d_i * Vg + Ve
      V_i <- d[i] * Vg + Ve  # n_traits x n_traits

      # Invert V_i (small matrix: n_traits x n_traits, typically 2-10)
      V_i_inv <- tryCatch(
        solve(V_i),
        error = function(e) solve(V_i + diag(1e-8, n_traits))
      )

      y_i <- Y_star[i, ]

      # Log-likelihood contribution: -0.5 * (log|V_i| + y_i' V_i^{-1} y_i)
      eig_Vi <- eigen(V_i, symmetric = TRUE, only.values = TRUE)$values
      Viy <- V_i_inv %*% y_i
      loglik <- loglik - 0.5 * (sum(log(pmax(eig_Vi, 1e-300))) +
                                   sum(y_i * Viy))

      # Pre-compute: V_i^{-1} y_i y_i' V_i^{-1}  (outer product in V^{-1} space)
      yyt_Vi <- tcrossprod(Viy)

      # ================================================================
      # Correct EM-REML updates (Meyer 1997; Lee & van der Werf 2006)
      #
      # For spectral model y*_i ~ N(0, d_i * Vg + Ve):
      #   E[a_i a_i'|y_i] = d_i*Vg - d_i^2*Vg*Vi^{-1}*Vg
      #                      + d_i^2*Vg*Vi^{-1}*y_i*y_i'*Vi^{-1}*Vg
      #
      #   Vg_new = (1/n) sum_i (1/d_i) * E[a_i a_i'|y_i]
      #          = (1/n) sum_i [Vg - d_i*Vg*Vi^{-1}*Vg
      #                         + d_i*Vg*Vi^{-1}*y_i*y_i'*Vi^{-1}*Vg]
      #
      #   E[e_i e_i'|y_i] = Ve - Ve*Vi^{-1}*Ve
      #                      + Ve*Vi^{-1}*y_i*y_i'*Vi^{-1}*Ve
      #
      #   Ve_new = (1/n) sum_i E[e_i e_i'|y_i]
      # ================================================================

      S_g <- Vg %*% V_i_inv  # Vg * V_i^{-1}

      Vg_acc <- Vg_acc + (Vg - d[i] * S_g %*% Vg +
                             d[i] * Vg %*% yyt_Vi %*% Vg)

      S_e <- Ve %*% V_i_inv  # Ve * V_i^{-1}

      Ve_acc <- Ve_acc + (Ve - S_e %*% Ve +
                             Ve %*% yyt_Vi %*% Ve)
    }

    Vg_em <- Vg_acc / n
    Ve_em <- Ve_acc / n

    # ---- Step-size damping to prevent oscillation ----
    Vg <- alpha * Vg_em + (1 - alpha) * Vg_old
    Ve <- alpha * Ve_em + (1 - alpha) * Ve_old

    # Ensure symmetry and positive definite
    Vg <- (Vg + t(Vg)) / 2
    Ve <- (Ve + t(Ve)) / 2
    Vg <- make_pd(Vg)
    Ve <- make_pd(Ve)

    # ---- Convergence check on relative log-likelihood change ----
    if (is.finite(loglik_old) && is.finite(loglik) && abs(loglik_old) > 1e-10) {
      rel_change <- abs((loglik - loglik_old) / loglik_old)
    } else {
      rel_change <- Inf
    }

    if (verbose && iter %% 10 == 0) {
      cat(sprintf("  Iter %d : loglik = %.4f, rel_change = %.2e\n",
                  iter, loglik, rel_change))
    }

    if (is.finite(rel_change) && rel_change < tol && iter > 5) {
      converged <- TRUE
      final_iter <- iter
      if (verbose) cat("  Converged at iteration", iter,
                       "(loglik =", round(loglik, 4), ")\n")
      break
    }
  }

  if (!converged && verbose) {
    cat("  WARNING: EM-REML did not converge after", nIters,
        "iterations (rel_change =", sprintf("%.2e", rel_change), ")\n")
  }

  return(list(
    Vg = Vg,
    Ve = Ve,
    converged = converged,
    iterations = final_iter,
    loglik = loglik
  ))
}

#' Make Matrix Positive Definite
#' @param M Square matrix
#' @param tol Minimum eigenvalue tolerance
#' @return Positive definite matrix
#' @keywords internal
make_pd <- function(M, tol = 1e-8) {
  # Use the main make_positive_definite function
  make_positive_definite(M, tol = tol)
}

#' Safe cov2cor
#' @param V Covariance matrix
#' @return Correlation matrix
#' @keywords internal
cov2cor_safe <- function(V) {
  d <- sqrt(diag(V))
  d[d == 0] <- 1
  V / outer(d, d)
}

#' Calculate Selection Index
#'
#' Computes a linear selection index from multi-trait BLUPs using user-specified
#' economic weights. The index ranks genotypes based on aggregate merit:
#' I = w1*GEBV1 + w2*GEBV2 + ... + wt*GEBVt
#'
#' @param mt_model Multi-trait model from fit_multi_trait()
#' @param weights Named vector of weights for each trait
#' @return Data frame with GID and Index columns, sorted by descending index value
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
    tr <- traits[i]
    if (tr %in% names(weights)) {
      index_values <- index_values + weights[tr] * blups[, tr]
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
    for (tr in names(mt_model$heritability)) {
      cat("  ", tr, ":", round(mt_model$heritability[tr], 3), "\n")
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
    for (tr in mt_model$traits) {
      cat("  ", tr, ":", round(mt_model$genetic_variance[tr, tr], 3), "\n")
    }
    cat("\n")
  }

  if (!is.null(mt_model$residual_variance)) {
    cat("Residual Variances (diagonal):\n")
    for (tr in mt_model$traits) {
      cat("  ", tr, ":", round(mt_model$residual_variance[tr, tr], 3), "\n")
    }
  }

  invisible(NULL)
}
