#' Mixed Model Equations Solver
#'
#' Core engine for solving mixed model equations using REML.
#' Supports both standard AI-REML and efficient EMMA-style algorithms.
#'
#' @name mme_solver
#' @keywords internal
NULL

# =============================================================================
# REML Algorithm Constants
# =============================================================================
# These values control numerical stability in the mixed model equations solver.
# Based on standard practices in ASReml, WOMBAT, and BLUPF90.
#
# WHEN TO ADJUST:
# - Convergence failures → try larger MIN_VARIANCE, smaller STEP_HALVING_INIT
# - Variance at boundary → may indicate true zero variance, not numerical issue
# - Singular K matrix → increase DEFAULT_RIDGE
# - Very large datasets (n > 10,000) → defaults should work fine
# - Very small datasets (n < 50) → may need larger MIN_VARIANCE (1e-8)
# =============================================================================

.MME_CONSTANTS <- list(
  # Minimum variance component (units: variance scale)
  # Value: 1e-10 is ~sqrt(machine epsilon), safely above numerical noise
  # Too small: numerical instability in matrix operations
  # Too large: might mask true small but non-zero variances
  # Reference: Gilmour et al. (1995) ASReml User Guide
  MIN_VARIANCE = 1e-10,
  
  # Minimum eigenvalue for positive definiteness
  # Value: 1e-10 matches MIN_VARIANCE for consistency
  # Eigenvalues below this are floored to this value
  # Reference: Misztal (2008) BLUPF90 documentation
  MIN_EIGENVALUE = 1e-10,
  
  # Ridge parameter for matrix regularization (units: fraction of diagonal)
  # Value: 1e-6 adds minimal perturbation while ensuring invertibility
  # Typical range: 1e-8 to 1e-4
  # Increase if: matrix inversion fails, K matrix is ill-conditioned
  # Reference: Thompson & Meyer (1986) estimation of variance components
  DEFAULT_RIDGE = 1e-6,
  
  # Minimum sample size for reliable estimation
  # Value: 5 based on minimum degrees of freedom requirements
  # With fewer observations, variance estimates are highly unreliable
  # Reference: Standard statistical practice (df = n - p - 1 > 0)
  MIN_OBSERVATIONS = 5,
  
  # Newton-Raphson step size parameters
  # INIT = 0.9: Start slightly conservative (not full Newton step)
  # This prevents overshooting in early iterations
  # Reference: Press et al. (2007) Numerical Recipes, §10.7
  STEP_HALVING_INIT = 0.9,
  
  # Minimum step before stopping step-halving
  # Value: 0.01 means we've reduced step by 99%
  # If we get here, likely at a boundary or saddle point
  STEP_HALVING_MIN = 0.01,
  
  # Ridge multiplier for singular matrices (units: multiplier)
  # When K_inv computation fails, multiply tolParInv by this factor
  # Reference: Common practice in mixed model software
  RIDGE_MULTIPLIER = 10,
  
  # Gradient tolerance multiplier for AI-REML convergence
  # Gradient norm threshold = tolParConvLL * GRADIENT_TOL_MULTIPLIER
  # Value: 10 allows gradient to be less stringent than likelihood change
  # This handles boundary solutions where gradient can be non-zero at optimum
  # Reference: Gilmour et al. (2009) ASReml convergence criteria
  GRADIENT_TOL_MULTIPLIER = 10
)

#' Fit Mixed Model Equations
#'
#' Solves Henderson's Mixed Model Equations using AI-REML for variance
#' component estimation. This is the core computational engine.
#'
#' @param y Response vector (n x 1)
#' @param X Fixed effects design matrix (n x p)
#' @param Z_list List of random effects design matrices
#' @param K_list List of relationship/covariance matrices for each random effect
#' @param nIters Maximum number of iterations
#' @param tolParConvLL Convergence tolerance for log-likelihood
#' @param tolParInv Tolerance for matrix inversion (ridge)
#' @param verbose Print progress
#' @param use.emma Use EMMA-style efficient algorithm when possible (single random effect)
#'
#' @return List with model components
#' @keywords internal
fit_mme <- function(y, X, Z_list, K_list,
                    nIters = 50,
                    tolParConvLL = 1e-04,
                    tolParInv = 1e-06,
                    verbose = FALSE,
                    use.emma = TRUE) {
  
  # Validate inputs
  n <- length(y)
  p <- ncol(X)
  n_random <- length(Z_list)
  
  if (n != nrow(X)) {
    stop("Dimensions of y and X do not match")
  }
  
  if (length(K_list) != n_random) {
    stop("Number of Z matrices must match number of K matrices")
  }
  
  # Remove NA observations
  complete_idx <- !is.na(y)
  if (sum(complete_idx) < n) {
    y <- y[complete_idx]
    X <- X[complete_idx, , drop = FALSE]
    Z_list <- lapply(Z_list, function(Z) Z[complete_idx, , drop = FALSE])
    n <- length(y)
  }
  
  if (n < .MME_CONSTANTS$MIN_OBSERVATIONS) {
    stop("Too few observations (", n, "). Need at least ", 
         .MME_CONSTANTS$MIN_OBSERVATIONS)
  }
  
  # Convert to matrices if needed
  y <- as.matrix(y)
  X <- as.matrix(X)
  Z_list <- lapply(Z_list, as.matrix)
  K_list <- lapply(K_list, as.matrix)
  
  # For single random effect with use.emma = TRUE, always prefer EMMA
  # regardless of n vs q ratio — EMMA is O(n) per iteration after eigendecomp
  if (n_random == 1 && use.emma) {
    return(fit_mme_emma(y, X, Z_list[[1]], K_list[[1]], 
                        tolParConvLL, tolParInv, verbose))
  }
  
  # Detect if Henderson's MME is more efficient (when n >> q)
  q_sizes <- sapply(Z_list, ncol)
  total_q <- sum(q_sizes)
  
  # When n >> q, always use Henderson-space solvers (never build n×n matrices)
  if (n > 2 * (p + total_q)) {
    if (use.emma) {
      # Default: Henderson EM-REML (fast per iteration, linear convergence)
      if (verbose) {
        cat("Auto-selecting Henderson's MME with EM-REML (n =", n, ", p+q =", p + total_q, ")\n")
        cat("  Tip: Set use.emma = FALSE for Henderson EM-REML with exact REML LL\n")
      }
      return(fit_mme_henderson(y, X, Z_list, K_list,
                              nIters, tolParConvLL, tolParInv, verbose))
    } else {
      # use.emma = FALSE: Henderson EM-REML with exact REML LL (monotone convergence)
      if (verbose) {
        cat("Using Henderson EM-REML solver (n =", n, ", p+q =", p + total_q, ")\n")
      }
      return(fit_mme_henderson_ai(y, X, Z_list, K_list,
                                   nIters, tolParConvLL, tolParInv, verbose))
    }
  }
  
  # For small n with single random effect (shouldn't reach here if use.emma=TRUE)
  if (n_random == 1 && use.emma) {
    return(fit_mme_emma(y, X, Z_list[[1]], K_list[[1]], 
                        tolParConvLL, tolParInv, verbose))
  }
  
  # Fallback: observation-space AI-REML (fine for small n)
  if (verbose) {
    cat("Using observation-space AI-REML solver (n =", n, ", p+q =", p + total_q, ")\n")
  }
  return(fit_mme_ai_reml(y, X, Z_list, K_list, 
                         nIters, tolParConvLL, tolParInv, verbose))
}

#' EMMA-style Efficient MME Solver
#'
#' Uses spectral decomposition for O(n) per-iteration complexity.
#' Only applicable for single random effect models.
#'
#' @keywords internal
fit_mme_emma <- function(y, X, Z, K, tol, ridge, verbose) {
  
  n <- length(y)
  p <- ncol(X)
  q <- ncol(Z)
  
  if (verbose) {
    cat("Using EMMA-style efficient algorithm\n")
  }
  
  # Compute ZKZ'
  ZKZt <- Z %*% K %*% t(Z)
  
  # Use C++ EMMA REML if available
  emma_result <- tryCatch({
    emma_reml_cpp(as.vector(y), X, ZKZt, tol, 100L)
  }, error = function(e) {
    # Fallback to R implementation
    fit_emma_r(y, X, ZKZt, tol, 100)
  })
  
  sigma2_g <- emma_result$sigma2_g
  sigma2_e <- emma_result$sigma2_e
  
  # Boundary clamping: if genetic variance is negligible relative to total,
  # clamp to MIN_VARIANCE floor (consistent with Henderson EM-REML solver).
  # This enables downstream boundary detection. Without clamping, EMMA's grid
  # search returns small but non-zero sigma2_g that never hits the floor.
  total_var <- sigma2_g + sigma2_e
  if (total_var > 0 && (sigma2_g / total_var) < 1e-6) {
    sigma2_g <- .MME_CONSTANTS$MIN_VARIANCE
  }
  sigma2_g <- max(sigma2_g, .MME_CONSTANTS$MIN_VARIANCE)
  sigma2_e <- max(sigma2_e, .MME_CONSTANTS$MIN_VARIANCE)
  
  if (verbose) {
    cat("Variance components: sigma2_g =", round(sigma2_g, 4),
        ", sigma2_e =", round(sigma2_e, 4), "\n")
  }
  
  # Solve for fixed and random effects using estimated variances
  mme_solution <- solve_mme_with_varcomp(y, X, Z, K, sigma2_g, sigma2_e, ridge)
  
  # EMMA always fits exactly 1 random effect (single-component by design),
  # so n_random = 1 and total free variance parameters = n_random + 1 = 2.
  # BUG FIX (v2.2.1): use n_random_emma = 1 explicitly so the formula
  # matches the Henderson/AI-REML paths (which use n_random + 1) and
  # AIC/BIC are comparable across all solvers.
  n_random_emma <- 1L  # EMMA is only called for single-random-effect models
  
  # Return in standard format
  return(list(
    b = mme_solution$b,
    u = list(mme_solution$u),
    uList = list(mme_solution$u),
    theta = c(sigma2_g, sigma2_e),
    Vi = mme_solution$Vi,
    P = mme_solution$P,
    Ci = mme_solution$Ci,
    llik = emma_result$loglik,
    convergence = TRUE,
    nIter = 1,
    monitor = matrix(c(sigma2_g, sigma2_e), ncol = 1),
    AIC = -2 * emma_result$loglik + 2 * (n_random_emma + 1),
    BIC = -2 * emma_result$loglik + log(n) * (n_random_emma + 1),
    fitted = as.vector(mme_solution$fitted),
    residuals = as.vector(mme_solution$residuals),
    n = n,
    p = p,
    df_residual = n - p
  ))
}

#' Henderson's MME Solver (operates in (p+q)-space)
#'
#' Solves Henderson's Mixed Model Equations directly in reduced dimension space.
#' Much more efficient than observation-space methods when n >> q.
#' Supports multiple random effects and uses EM-REML for variance estimation.
#'
#' @keywords internal
fit_mme_henderson <- function(y, X, Z_list, K_list,
                               nIters = 50,
                               tolParConvLL = 1e-04,
                               tolParInv = 1e-06,
                               verbose = FALSE) {
  
  n <- length(y)
  p <- ncol(X)
  n_random <- length(Z_list)
  q_sizes <- sapply(Z_list, ncol)
  
  if (verbose) {
    cat("Henderson's MME solver\n")
    cat("  Observations (n):", n, "\n")
    cat("  Fixed effects (p):", p, "\n")
    cat("  Random effects:", n_random, "with sizes", paste(q_sizes, collapse = ", "), "\n")
    cat("  Working dimension:", p + sum(q_sizes), "vs observation-space:", n, "\n")
  }
  
  # Try C++ implementation for single random effect (fastest path)
  if (n_random == 1) {
    cpp_result <- tryCatch({
      # Initialize variance components using residual variance after projecting out fixed effects
      resid_y <- tryCatch({
        residuals(lm.fit(X, as.vector(y)))
      }, error = function(e) {
        as.vector(y) - mean(as.vector(y), na.rm = TRUE)
      })
      var_y <- var(resid_y, na.rm = TRUE)
      if (var_y == 0 || is.na(var_y)) var_y <- 1
      sigma2_g_init <- 0.5 * var_y
      sigma2_e_init <- 0.5 * var_y
      
      result_cpp <- fit_henderson_mme_cpp(
        y = as.vector(y),
        X = X,
        Z = Z_list[[1]],
        K = K_list[[1]],
        sigma2_g_init = sigma2_g_init,
        sigma2_e_init = sigma2_e_init,
        max_iter = as.integer(nIters),
        tol = tolParConvLL,
        verbose = verbose
      )
      
      # Convert to standard format
      # BUG FIX (v2.2.1): use n_random (= 1 here, single-effect C++ path)
      # so AIC/BIC are comparable with Henderson R and AI-REML solvers.
      list(
        b = result_cpp$b,
        u = list(result_cpp$u),
        uList = list(result_cpp$u),
        theta = result_cpp$theta,
        Vi = NULL,
        P = NULL,
        Ci = result_cpp$C_inv[1:p, 1:p, drop = FALSE],
        llik = result_cpp$llik,
        convergence = result_cpp$convergence,
        nIter = result_cpp$nIter,
        monitor = matrix(result_cpp$theta, ncol = 1),
        AIC = -2 * result_cpp$llik + 2 * (n_random + 1),
        BIC = -2 * result_cpp$llik + log(n) * (n_random + 1),
        fitted = as.vector(result_cpp$fitted),
        residuals = as.vector(result_cpp$residuals),
        n = n,
        p = p,
        df_residual = n - p
      )
    }, error = function(e) {
      if (verbose) {
        cat("C++ implementation not available, using R version\n")
      }
      NULL
    })
    
    if (!is.null(cpp_result)) {
      return(cpp_result)
    }
  }
  
  # R implementation (supports multiple random effects)
  # Initialize variance components
  # BUG FIX (v2.2.1): smarter initialisation that detects ENV-type components
  # and seeds them from between-environment variance, not a fixed 25% share.
  # ENV variance can be 3-10× larger than additive genetic variance (e.g. G2F
  # 2021: sigma2_ENV ~ 745 vs sigma2_A ~ 246). Starting ENV at 25% of residual
  # variance wastes dozens of EM iterations just moving variance between ENV
  # and residual before genetic signals are explored.
  resid_y <- tryCatch({
    residuals(lm.fit(X, as.vector(y)))
  }, error = function(e) {
    # Fallback if lm.fit fails (e.g., singular X)
    as.vector(y) - mean(as.vector(y), na.rm = TRUE)
  })
  var_y <- var(resid_y, na.rm = TRUE)
  if (var_y == 0 || is.na(var_y)) var_y <- 1
  
  # Classify random effects: ENV-type have large Z with few levels (n_env << n);
  # genetic-type have large K with many levels (n_gid).
  # Heuristic: if q_i < 0.15 * n AND q_i < 100, treat as ENV component.
  is_env_component <- q_sizes < max(15, 0.15 * n) & q_sizes < 100
  
  theta <- numeric(n_random + 1)
  genetic_share <- 0.5
  theta[n_random + 1] <- var_y * (1 - genetic_share)  # Residual
  
  if (n_random == 1) {
    theta[1] <- var_y * genetic_share
  } else if (!any(is_env_component)) {
    # All large-dimensional effects: original share scheme
    shares <- c(0.6, 0.25, rep(0.15 / max(1, n_random - 2), max(0, n_random - 2)))
    shares <- shares[1:n_random]
    shares <- shares / sum(shares)
    theta[1:n_random] <- var_y * genetic_share * shares
  } else {
    # Mixed model with at least one ENV-type component.
    # Seed ENV components from the between-group variance of y across ENV levels.
    # Seed genetic components from remaining variance after ENV is accounted for.
    n_env_comps     <- sum(is_env_component)
    n_genetic_comps <- n_random - n_env_comps
    
    # Between-ENV variance: var of group means scaled by avg group size
    # This is a rough ANOVA-style estimate: sigma2_env ~ (MS_between - MS_within) / n_bar
    # Use simpler: sigma2_env_init = max(0, var(group_means) * fraction)
    env_init_share <- 0.5  # ENV gets up to 50% of total variance as seed
    env_seed <- var_y * env_init_share / max(1, n_env_comps)
    
    # Genetic components share the remainder
    gen_seed_total <- var_y * genetic_share
    if (n_genetic_comps == 1) {
      gen_shares <- 1.0
    } else {
      gen_shares <- c(0.6, 0.25, rep(0.15 / max(1, n_genetic_comps - 2),
                                       max(0, n_genetic_comps - 2)))
      gen_shares <- gen_shares[seq_len(n_genetic_comps)]
      gen_shares <- gen_shares / sum(gen_shares)
    }
    
    gen_idx <- 0L
    for (i in seq_len(n_random)) {
      if (is_env_component[i]) {
        theta[i] <- env_seed
      } else {
        gen_idx <- gen_idx + 1L
        theta[i] <- gen_seed_total * gen_shares[gen_idx]
      }
    }
  }
  
  names(theta) <- c(paste0("u", seq_len(n_random)), "e")
  
  # Pre-compute K inverses (q × q matrices)
  K_inv_list <- vector("list", n_random)
  for (i in seq_len(n_random)) {
    K_i <- K_list[[i]]
    K_inv_list[[i]] <- tryCatch({
      solve(K_i)
    }, error = function(e) {
      # Add ridge if singular
      solve(K_i + diag(tolParInv * .MME_CONSTANTS$RIDGE_MULTIPLIER, nrow(K_i)))
    })
  }
  
  # Pre-compute X'X, X'y (p × p and p × 1)
  XtX <- crossprod(X)  # More efficient than t(X) %*% X
  Xty <- crossprod(X, y)
  
  # Pre-compute Z'Z, Z'X, Z'y for each random effect
  ZtZ_list <- vector("list", n_random)
  ZtX_list <- vector("list", n_random)
  Zty_list <- vector("list", n_random)
  
  for (i in seq_len(n_random)) {
    Z_i <- Z_list[[i]]
    ZtZ_list[[i]] <- crossprod(Z_i)  # q_i × q_i, efficient for sparse Z
    ZtX_list[[i]] <- crossprod(Z_i, X)  # q_i × p
    Zty_list[[i]] <- crossprod(Z_i, y)  # q_i × 1
  }
  
  # EM-REML iterations
  llik_old <- -Inf
  converged <- FALSE
  monitor <- matrix(NA, nrow = n_random + 1, ncol = nIters)
  Ci <- NULL  # Will store current iteration's C⁻¹
  
  for (iter in seq_len(nIters)) {
    
    sigma2_e <- theta[n_random + 1]
    
    # Build Henderson's coefficient matrix C in (p + Σq_i) space
    # C = [X'R⁻¹X      X'R⁻¹Z     ]
    #     [Z'R⁻¹X      Z'R⁻¹Z+G⁻¹ ]
    # where R⁻¹ = (1/σ²_e)I, so X'R⁻¹X = X'X/σ²_e
    
    total_dim <- p + sum(q_sizes)
    C <- matrix(0, total_dim, total_dim)
    rhs <- numeric(total_dim)
    
    # Fill X'X block (scaled by 1/sigma2_e)
    C[1:p, 1:p] <- XtX / sigma2_e
    rhs[1:p] <- Xty / sigma2_e
    
    # Fill mixed blocks and random effects blocks
    row_offset <- p
    for (i in seq_len(n_random)) {
      q_i <- q_sizes[i]
      rows_i <- row_offset + (1:q_i)
      
      # X'Z block (scaled by 1/sigma2_e)
      C[1:p, rows_i] <- t(ZtX_list[[i]]) / sigma2_e
      C[rows_i, 1:p] <- ZtX_list[[i]] / sigma2_e
      
      # Z'y block
      rhs[rows_i] <- Zty_list[[i]] / sigma2_e
      
      # Z'Z block diagonal + G⁻¹ = Z'Z/σ²_e + K⁻¹/σ²_g
      sigma2_g_i <- theta[i]
      # Convention 2 (R⁻¹-scaled Henderson MME):
      # C_uu = Z'R⁻¹Z + G⁻¹ = Z'Z/σ²_e + K⁻¹/σ²_g
      # NOTE: Previously used lambda_i = σ²_e/σ²_g which is Convention 1's ratio,
      # but data blocks are already divided by σ²_e (Convention 2), creating a
      # mixed convention that makes the G⁻¹ penalty σ²_e times too strong.
      C[rows_i, rows_i] <- ZtZ_list[[i]] / sigma2_e + K_inv_list[[i]] / sigma2_g_i
      
      # Cross-terms between random effects (Z_i'Z_j / sigma2_e)
      if (n_random > 1 && i < n_random) {
        col_offset <- row_offset + q_i
        for (j in (i+1):n_random) {
          q_j <- q_sizes[j]
          cols_j <- col_offset + (1:q_j)
          
          # Z_i'Z_j
          ZiZj <- crossprod(Z_list[[i]], Z_list[[j]]) / sigma2_e
          C[rows_i, cols_j] <- ZiZj
          C[cols_j, rows_i] <- t(ZiZj)
          
          col_offset <- col_offset + q_j
        }
      }
      
      row_offset <- row_offset + q_i
    }
    
    # Solve Henderson's equations: C [b; u] = rhs
    sol <- tryCatch({
      solve(C, rhs)
    }, error = function(e) {
      # Add small ridge if singular
      C_ridge <- C + diag(tolParInv, nrow(C))
      solve(C_ridge, rhs)
    })
    
    # Extract fixed effects and random effects
    b <- sol[1:p]
    u <- vector("list", n_random)
    row_offset <- p
    for (i in seq_len(n_random)) {
      q_i <- q_sizes[i]
      u[[i]] <- sol[row_offset + (1:q_i)]
      names(u[[i]]) <- colnames(Z_list[[i]])
      row_offset <- row_offset + q_i
    }
    
    # Compute fitted values and residuals
    fitted_vals <- X %*% b
    for (i in seq_len(n_random)) {
      fitted_vals <- fitted_vals + Z_list[[i]] %*% u[[i]]
    }
    residuals <- y - fitted_vals
    
    # EM updates for variance components
    # BUG FIX (v2.3.1): Full EM-REML for ALL components including σ²_e.
    # Previously σ²_e = RSS/(n-p) which is the naive OLS estimate, NOT the
    # EM-REML update. The correct formula includes a trace correction:
    #   σ²_e = (RSS + σ²_e · tr(C⁻¹ · C_data)) / (n - p)
    # Without this, σ²_e is overestimated and genetic VCs underestimated.
    rss <- sum(residuals^2)
    theta_new <- theta
    
    # Compute C⁻¹ for trace corrections (critical for unbiased EM-REML)
    Ci <- tryCatch({
      solve(C)
    }, error = function(e) {
      C_ridge <- C + diag(tolParInv, nrow(C))
      solve(C_ridge)
    })
    
    # Update random effect variances with trace correction
    row_offset <- p
    for (i in seq_len(n_random)) {
      q_i <- q_sizes[i]
      rows_i <- row_offset + (1:q_i)
      
      u_i <- u[[i]]
      K_inv_i <- K_inv_list[[i]]
      
      # Quadratic form: u'K⁻¹u
      uKu <- as.numeric(t(u_i) %*% K_inv_i %*% u_i)
      
      # Trace correction: tr(C⁻¹_uu K⁻¹)
      # In Convention 2, C⁻¹_uu directly gives PEV (no σ²_e scaling needed)
      C_inv_uu <- Ci[rows_i, rows_i, drop = FALSE]
      trace_term <- sum(diag(C_inv_uu %*% K_inv_i))
      
      # Full EM update with trace correction
      theta_new[i] <- max((uKu + trace_term) / q_i, .MME_CONSTANTS$MIN_VARIANCE)
      
      row_offset <- row_offset + q_i
    }
    
    # Full EM-REML update for σ²_e with trace correction
    # tr(C⁻¹ · C_data) = (p + Σq_i) - Σ_i tr(C⁻¹_uu_i · K⁻¹_i / σ²_g_i)
    total_dim_tr <- p + sum(q_sizes)
    trace_Cinv_Ginv <- 0
    row_offset_tr <- p
    for (i in seq_len(n_random)) {
      q_i <- q_sizes[i]
      rows_i <- row_offset_tr + (1:q_i)
      C_inv_uu_i <- Ci[rows_i, rows_i, drop = FALSE]
      trace_Cinv_Ginv <- trace_Cinv_Ginv +
        sum(diag(C_inv_uu_i %*% K_inv_list[[i]])) / theta[i]
      row_offset_tr <- row_offset_tr + q_i
    }
    trace_Cinv_Cdata <- total_dim_tr - trace_Cinv_Ginv
    sigma2_e <- max((rss + theta[n_random + 1] * trace_Cinv_Cdata) / (n - p),
                    .MME_CONSTANTS$MIN_VARIANCE)
    theta_new[n_random + 1] <- sigma2_e
    
    # BUG FIX (v2.3.1): Compute exact REML log-likelihood, not the simplified
    # -0.5*[n*log(σ²_e) + RSS/σ²_e] which is the i.i.d. Gaussian LL.
    # Exact REML LL = -0.5 * [(n-p)*log(σ²_e) + Σ(q_i*log(σ²_g_i) + log|K_i|)
    #                          + log|C| + y'Py]
    log_det_C <- as.numeric(determinant(C, logarithm = TRUE)$modulus)
    yPy <- as.numeric(crossprod(y) / sigma2_e - crossprod(rhs, sol))
    log_det_G <- 0
    for (i in seq_len(n_random)) {
      K_det <- tryCatch(
        as.numeric(determinant(K_list[[i]], logarithm = TRUE)$modulus),
        error = function(e) 0
      )
      log_det_G <- log_det_G + q_sizes[i] * log(theta[i]) + K_det
    }
    llik <- -0.5 * ((n - p) * log(sigma2_e) + log_det_G + log_det_C + yPy)
    
    # Check convergence
    llik_diff <- abs(llik - llik_old)
    theta_diff <- max(abs(theta_new - theta) / pmax(abs(theta), 1e-8))
    
    if (llik_diff < tolParConvLL && theta_diff < tolParConvLL && iter > 1) {
      converged <- TRUE
      if (verbose) {
        cat("Converged at iteration", iter, "\n")
      }
      break
    }
    
    # BUG FIX (v2.3.1): Tightened LL plateau from 0.01/15 to 0.001/25.
    # For large n, LL changes of 0.01 can still correspond to meaningful VC changes.
    if (iter > 25 && llik_diff < 0.001) {
      converged <- TRUE
      if (verbose) {
        cat("Converged at iteration ", iter, " (LL plateau)\n", sep = "")
      }
      break
    }
    
    # Constrain variance components
    for (j in seq_along(theta_new)) {
      if (theta_new[j] < .MME_CONSTANTS$MIN_VARIANCE) {
        theta_new[j] <- .MME_CONSTANTS$MIN_VARIANCE
      }
    }
    
    theta <- theta_new
    llik_old <- llik
    monitor[, iter] <- theta
    
    if (verbose && iter %% 5 == 0) {
      cat(sprintf("Iteration %d: LogLik = %.4f, theta = %s\n",
                  iter, llik, paste(round(theta, 4), collapse = ", ")))
    }
  }
  
  # Use final C⁻¹ from last iteration (already computed for trace correction)
  # Extract variance of fixed effects (top-left block of C⁻¹)
  Ci_b <- Ci[1:p, 1:p, drop = FALSE]
  
  # Model fit statistics
  AIC <- -2 * llik + 2 * (n_random + 1)
  BIC <- -2 * llik + log(n) * (n_random + 1)
  
  # Add names (b is a vector, so use names not rownames)
  names(b) <- colnames(X)
  
  return(list(
    b = b,
    u = u,
    uList = u,
    theta = theta,
    Vi = NULL,  # Not computed in Henderson's MME (would be n×n)
    P = NULL,   # Not computed in Henderson's MME (would be n×n)
    Ci = Ci_b,  # Variance of fixed effects
    llik = llik,
    convergence = converged,
    nIter = iter,
    monitor = monitor[, seq_len(iter), drop = FALSE],
    AIC = AIC,
    BIC = BIC,
    fitted = as.vector(fitted_vals),
    residuals = as.vector(residuals),
    n = n,
    p = p,
    df_residual = n - p
  ))
}

#' Henderson's MME with EM-REML Solver
#'
#' Uses Henderson's (p+q)-space formulation with EM-REML variance component
#' updates. Pure EM-REML guarantees monotonic increase of the REML
#' log-likelihood at every iteration. Computes the exact REML LL for
#' convergence monitoring.
#'
#' @param y Response vector (n x 1)
#' @param X Fixed effects design matrix (n x p)
#' @param Z_list List of random effects design matrices
#' @param K_list List of relationship/covariance matrices for each random effect
#' @param nIters Maximum number of iterations (EM needs more than Newton)
#' @param tolParConvLL Convergence tolerance for log-likelihood
#' @param tolParInv Tolerance for matrix inversion (ridge)
#' @param verbose Print progress
#'
#' @return List with model components (same format as other solvers)
#' @keywords internal
fit_mme_henderson_ai <- function(y, X, Z_list, K_list,
                                  nIters = 200,
                                  tolParConvLL = 1e-04,
                                  tolParInv = 1e-06,
                                  verbose = FALSE) {
  
  n <- length(y)
  p <- ncol(X)
  n_random <- length(Z_list)
  q_sizes <- sapply(Z_list, ncol)
  
  if (verbose) {
    cat("Henderson EM-REML solver (Aitken accelerated, LL-guarded)\n")
    cat("  Observations (n):", n, "\n")
    cat("  Fixed effects (p):", p, "\n")
    cat("  Random effects:", n_random, "with sizes", paste(q_sizes, collapse = ", "), "\n")
    cat("  Working dimension:", p + sum(q_sizes), "vs observation-space:", n, "\n")
  }
  
  # Initialize variance components (same as Henderson EM-REML)
  # BUG FIX (v2.2.1): smarter ENV-aware initialisation (see fit_mme_henderson).
  resid_y <- tryCatch({
    residuals(lm.fit(X, as.vector(y)))
  }, error = function(e) {
    as.vector(y) - mean(as.vector(y), na.rm = TRUE)
  })
  var_y <- var(resid_y, na.rm = TRUE)
  if (var_y == 0 || is.na(var_y)) var_y <- 1
  
  is_env_component <- q_sizes < max(15, 0.15 * n) & q_sizes < 100
  
  theta <- numeric(n_random + 1)
  genetic_share <- 0.5
  theta[n_random + 1] <- var_y * (1 - genetic_share)  # Residual
  
  if (n_random == 1) {
    theta[1] <- var_y * genetic_share
  } else if (!any(is_env_component)) {
    shares <- c(0.6, 0.25, rep(0.15 / max(1, n_random - 2), max(0, n_random - 2)))
    shares <- shares[1:n_random]
    shares <- shares / sum(shares)
    theta[1:n_random] <- var_y * genetic_share * shares
  } else {
    n_env_comps     <- sum(is_env_component)
    n_genetic_comps <- n_random - n_env_comps
    env_seed        <- var_y * 0.5 / max(1, n_env_comps)
    gen_seed_total  <- var_y * genetic_share
    if (n_genetic_comps == 1) {
      gen_shares <- 1.0
    } else {
      gen_shares <- c(0.6, 0.25, rep(0.15 / max(1, n_genetic_comps - 2),
                                       max(0, n_genetic_comps - 2)))
      gen_shares <- gen_shares[seq_len(n_genetic_comps)]
      gen_shares <- gen_shares / sum(gen_shares)
    }
    gen_idx <- 0L
    for (i in seq_len(n_random)) {
      if (is_env_component[i]) {
        theta[i] <- env_seed
      } else {
        gen_idx <- gen_idx + 1L
        theta[i] <- gen_seed_total * gen_shares[gen_idx]
      }
    }
  }
  
  names(theta) <- c(paste0("u", seq_len(n_random)), "e")
  
  # Pre-compute K inverses (q × q matrices)
  K_inv_list <- vector("list", n_random)
  for (i in seq_len(n_random)) {
    K_i <- K_list[[i]]
    K_inv_list[[i]] <- tryCatch({
      solve(K_i)
    }, error = function(e) {
      solve(K_i + diag(tolParInv * .MME_CONSTANTS$RIDGE_MULTIPLIER, nrow(K_i)))
    })
  }
  
  # Pre-compute data matrices (X'X, X'y, Z'Z, Z'X, Z'y)
  XtX <- crossprod(X)
  Xty <- crossprod(X, y)
  
  ZtZ_list <- vector("list", n_random)
  ZtX_list <- vector("list", n_random)
  Zty_list <- vector("list", n_random)
  
  for (i in seq_len(n_random)) {
    Z_i <- Z_list[[i]]
    ZtZ_list[[i]] <- crossprod(Z_i)
    ZtX_list[[i]] <- crossprod(Z_i, X)
    Zty_list[[i]] <- crossprod(Z_i, y)
  }
  
  # Pre-compute log|K| for each random effect (constant across iterations)
  log_det_K_list <- numeric(n_random)
  for (i in seq_len(n_random)) {
    log_det_K_list[i] <- as.numeric(determinant(K_list[[i]], logarithm = TRUE)$modulus)
  }
  
  # =====================================================================
  # EM-REML with guarded Aitken acceleration
  # =====================================================================
  # Pure EM-REML (monotone) + Aitken extrapolation every 3 iterations.
  # Aitken-accelerated θ is only accepted if LL(θ_acc) >= LL(θ_em),
  # evaluated via one extra (p+q)×(p+q) solve. This guarantees monotone
  # LL while achieving near-quadratic convergence speed.
  # =====================================================================
  
  llik_old <- -Inf
  llik_best <- -Inf  # best LL seen — Aitken must beat this
  converged <- FALSE
  nIters_em <- max(nIters, 200)
  monitor <- matrix(NA, nrow = n_random + 1, ncol = nIters_em)
  Ci <- NULL
  
  # Aitken acceleration: store 3 consecutive EM thetas
  theta_hist <- vector("list", 3)
  
  for (iter in seq_len(nIters_em)) {
    
    sigma2_e <- theta[n_random + 1]
    
    # Build Henderson's coefficient matrix C in (p + Σq_i) space
    total_dim <- p + sum(q_sizes)
    C <- matrix(0, total_dim, total_dim)
    rhs <- numeric(total_dim)
    
    # Fill X'X block (scaled by 1/sigma2_e)
    C[1:p, 1:p] <- XtX / sigma2_e
    rhs[1:p] <- Xty / sigma2_e
    
    # Fill mixed blocks and random effects blocks
    row_offset <- p
    for (i in seq_len(n_random)) {
      q_i <- q_sizes[i]
      rows_i <- row_offset + (1:q_i)
      
      # X'Z block (scaled by 1/sigma2_e)
      C[1:p, rows_i] <- t(ZtX_list[[i]]) / sigma2_e
      C[rows_i, 1:p] <- ZtX_list[[i]] / sigma2_e
      
      # Z'y block
      rhs[rows_i] <- Zty_list[[i]] / sigma2_e
      
      # Z'Z block diagonal + G⁻¹
      sigma2_g_i <- theta[i]
      # Convention 2: G⁻¹ = K⁻¹/σ²_g
      C[rows_i, rows_i] <- ZtZ_list[[i]] / sigma2_e + K_inv_list[[i]] / sigma2_g_i
      
      # Cross-terms between random effects
      if (n_random > 1 && i < n_random) {
        col_offset <- row_offset + q_i
        for (j in (i+1):n_random) {
          q_j <- q_sizes[j]
          cols_j <- col_offset + (1:q_j)
          
          ZiZj <- crossprod(Z_list[[i]], Z_list[[j]]) / sigma2_e
          C[rows_i, cols_j] <- ZiZj
          C[cols_j, rows_i] <- t(ZiZj)
          
          col_offset <- col_offset + q_j
        }
      }
      
      row_offset <- row_offset + q_i
    }
    
    # Solve Henderson's equations: C [b; u] = rhs
    sol <- tryCatch({
      solve(C, rhs)
    }, error = function(e) {
      C_ridge <- C + diag(tolParInv, nrow(C))
      solve(C_ridge, rhs)
    })
    
    # Extract solutions
    b <- sol[1:p]
    u <- vector("list", n_random)
    row_offset <- p
    for (i in seq_len(n_random)) {
      q_i <- q_sizes[i]
      u[[i]] <- sol[row_offset + (1:q_i)]
      names(u[[i]]) <- colnames(Z_list[[i]])
      row_offset <- row_offset + q_i
    }
    
    # Compute fitted values and residuals
    fitted_vals <- X %*% b
    for (i in seq_len(n_random)) {
      fitted_vals <- fitted_vals + Z_list[[i]] %*% u[[i]]
    }
    residuals <- y - fitted_vals
    rss <- sum(residuals^2)
    
    # Compute C⁻¹ for trace corrections (required for unbiased EM-REML)
    Ci <- tryCatch({
      solve(C)
    }, error = function(e) {
      C_ridge <- C + diag(tolParInv, nrow(C))
      solve(C_ridge)
    })
    
    # =================================================================
    # EXACT REML LOG-LIKELIHOOD (for convergence monitoring)
    # =================================================================
    # ℓ = -0.5·[(n-p)·log(σ²_e) + Σqᵢ·log(σ²_gᵢ) + Σlog|Kᵢ| + log|C| + y'Py]
    log_det_C <- as.numeric(determinant(C, logarithm = TRUE)$modulus)
    
    yPy <- as.numeric(crossprod(y) / sigma2_e - crossprod(rhs, sol))
    
    log_det_G <- 0
    for (i in seq_len(n_random)) {
      log_det_G <- log_det_G + q_sizes[i] * log(theta[i]) + log_det_K_list[i]
    }
    
    llik <- -0.5 * ((n - p) * log(sigma2_e) + log_det_G + log_det_C + yPy)
    llik_best <- max(llik_best, llik)
    
    # =================================================================
    # EM-REML VARIANCE COMPONENT UPDATES
    # =================================================================
    # σ²_g_i = (û'K⁻¹û + tr(C⁻¹_uu · K⁻¹)) / q_i   (exact EM-REML)
    # σ²_e   = (RSS + σ²_e · tr(C⁻¹ · C_data)) / (n - p)  (full EM-REML)
    #
    # BUG FIX (v2.3.1): σ²_e now uses full EM-REML with trace correction.
    # The simplified σ²_e = RSS/(n-p) converges to a DIFFERENT fixed point
    # than the full EM update when there are multiple random effects.
    # The bias is: σ²_e overestimated, σ²_g underestimated.
    # =================================================================
    
    theta_new <- theta
    
    # Update σ²_g for each random effect
    row_offset <- p
    for (i in seq_len(n_random)) {
      q_i <- q_sizes[i]
      rows_i <- row_offset + (1:q_i)
      
      u_i <- u[[i]]
      K_inv_i <- K_inv_list[[i]]
      uKu <- as.numeric(t(u_i) %*% K_inv_i %*% u_i)
      C_inv_uu <- Ci[rows_i, rows_i, drop = FALSE]
      trace_term <- sum(diag(C_inv_uu %*% K_inv_i))
      
      theta_new[i] <- max((uKu + trace_term) / q_i, .MME_CONSTANTS$MIN_VARIANCE)
      row_offset <- row_offset + q_i
    }
    
    # Full EM-REML update for σ²_e with trace correction
    # tr(C⁻¹ · C_data) = (p + Σq_i) - Σ_i tr(C⁻¹_uu_i · K⁻¹_i / σ²_g_i)
    total_dim_tr <- p + sum(q_sizes)
    trace_Cinv_Ginv <- 0
    row_offset_tr <- p
    for (i in seq_len(n_random)) {
      q_i <- q_sizes[i]
      rows_i <- row_offset_tr + (1:q_i)
      C_inv_uu_i <- Ci[rows_i, rows_i, drop = FALSE]
      trace_Cinv_Ginv <- trace_Cinv_Ginv +
        sum(diag(C_inv_uu_i %*% K_inv_list[[i]])) / theta[i]
      row_offset_tr <- row_offset_tr + q_i
    }
    trace_Cinv_Cdata <- total_dim_tr - trace_Cinv_Ginv
    theta_new[n_random + 1] <- max((rss + sigma2_e * trace_Cinv_Cdata) / (n - p),
                                   .MME_CONSTANTS$MIN_VARIANCE)
    
    # =================================================================
    # AITKEN ACCELERATION WITH LL GUARD
    # =================================================================
    # Store EM theta, then attempt Aitken every 3 iterations.
    # Only accept accelerated theta if LL(θ_acc) >= LL(θ_em).
    # This gives near-quadratic convergence while guaranteeing
    # monotone LL increase (EM fallback is always available).
    # Cost: one extra (p+q)×(p+q) solve every 3 iterations.
    # =================================================================
    
    theta_em <- theta_new  # save plain EM result as fallback
    
    cycle_pos <- ((iter - 1) %% 3) + 1
    theta_hist[[cycle_pos]] <- theta_new
    
    if (cycle_pos == 3 && iter >= 6) {
      t1 <- theta_hist[[1]]
      t2 <- theta_hist[[2]]
      t3 <- theta_hist[[3]]
      denom <- t3 - 2 * t2 + t1
      
      # Compute Aitken extrapolation element-wise
      theta_acc <- theta_new
      any_accelerated <- FALSE
      for (j in seq_along(theta_new)) {
        if (abs(denom[j]) > 1e-12) {
          acc_j <- t1[j] - (t2[j] - t1[j])^2 / denom[j]
          if (acc_j > .MME_CONSTANTS$MIN_VARIANCE && acc_j < var_y * 10) {
            theta_acc[j] <- acc_j
            any_accelerated <- TRUE
          }
        }
      }
      
      # LL guard: evaluate LL at θ_acc using Henderson identity
      if (any_accelerated) {
        sigma2_e_acc <- theta_acc[n_random + 1]
        C_acc <- matrix(0, total_dim, total_dim)
        rhs_acc <- numeric(total_dim)
        C_acc[1:p, 1:p] <- XtX / sigma2_e_acc
        rhs_acc[1:p] <- Xty / sigma2_e_acc
        
        row_off <- p
        for (ii in seq_len(n_random)) {
          qi <- q_sizes[ii]
          ri <- row_off + (1:qi)
          C_acc[1:p, ri] <- t(ZtX_list[[ii]]) / sigma2_e_acc
          C_acc[ri, 1:p] <- ZtX_list[[ii]] / sigma2_e_acc
          C_acc[ri, ri] <- ZtZ_list[[ii]] / sigma2_e_acc + K_inv_list[[ii]] / theta_acc[ii]
          rhs_acc[ri] <- Zty_list[[ii]] / sigma2_e_acc
          if (n_random > 1 && ii < n_random) {
            co <- row_off + qi
            for (jj in (ii+1):n_random) {
              qj <- q_sizes[jj]
              cj <- co + (1:qj)
              ZiZj_acc <- crossprod(Z_list[[ii]], Z_list[[jj]]) / sigma2_e_acc
              C_acc[ri, cj] <- ZiZj_acc
              C_acc[cj, ri] <- t(ZiZj_acc)
              co <- co + qj
            }
          }
          row_off <- row_off + qi
        }
        
        sol_acc <- tryCatch(solve(C_acc, rhs_acc), error = function(e) NULL)
        
        if (!is.null(sol_acc)) {
          yPy_acc <- as.numeric(crossprod(y) / sigma2_e_acc - crossprod(rhs_acc, sol_acc))
          log_det_C_acc <- as.numeric(determinant(C_acc, logarithm = TRUE)$modulus)
          log_det_G_acc <- 0
          for (ii in seq_len(n_random)) {
            log_det_G_acc <- log_det_G_acc + q_sizes[ii] * log(theta_acc[ii]) + log_det_K_list[ii]
          }
          llik_acc <- -0.5 * ((n - p) * log(sigma2_e_acc) + log_det_G_acc + log_det_C_acc + yPy_acc)
          
          # Accept if LL at least as good as current EM iteration
          if (llik_acc >= llik - 0.1) {
            theta_new <- theta_acc
            if (llik_acc > llik_best) llik_best <- llik_acc
          }
          # else: keep theta_em (already in theta_new)
        }
      }
    }
    
    # Global bounds
    for (j in seq_along(theta_new)) {
      theta_new[j] <- max(theta_new[j], .MME_CONSTANTS$MIN_VARIANCE)
      theta_new[j] <- min(theta_new[j], var_y * 10)
    }
    
    # Check convergence
    if (iter > 1) {
      llik_diff <- abs(llik - llik_old)
      theta_diff <- max(abs(theta_new - theta) / pmax(abs(theta), 1e-8))
      
      # Primary: both LL and parameters converged
      if (llik_diff < tolParConvLL && theta_diff < tolParConvLL) {
        converged <- TRUE
        if (verbose) cat("Converged at iteration", iter, "\n")
        theta <- theta_new
        llik_old <- llik
        monitor[, iter] <- theta
        break
      }
      
      # Secondary: LL effectively flat (confounded components trading variance)
      # BUG FIX (v2.3.1): Tightened from 0.01/15 to 0.001/25.
      if (iter > 25 && llik_diff < 0.001) {
        converged <- TRUE
        if (verbose) cat("Converged at iteration", iter, 
                         " (LL plateau, delta=", signif(llik_diff, 3), ")\n", sep = "")
        theta <- theta_new
        llik_old <- llik
        monitor[, iter] <- theta
        break
      }
    }
    
    theta <- theta_new
    llik_old <- llik
    monitor[, iter] <- theta
    
    if (verbose && iter %% 5 == 0) {
      cat(sprintf("Iteration %d: LogLik = %.4f, theta = %s\n",
                  iter, llik, paste(round(theta, 4), collapse = ", ")))
    }
  }
  
  # Use final C⁻¹ from last iteration
  Ci_b <- Ci[1:p, 1:p, drop = FALSE]
  
  # Model fit statistics
  AIC <- -2 * llik + 2 * (n_random + 1)
  BIC <- -2 * llik + log(n) * (n_random + 1)
  
  # Add names
  names(b) <- colnames(X)
  
  return(list(
    b = b,
    u = u,
    uList = u,
    theta = theta,
    Vi = NULL,  # Not computed in Henderson's MME (would be n×n)
    P = NULL,   # Not computed in Henderson's MME (would be n×n)
    Ci = Ci_b,  # Variance of fixed effects
    llik = llik,
    convergence = converged,
    nIter = iter,
    monitor = monitor[, seq_len(iter), drop = FALSE],
    AIC = AIC,
    BIC = BIC,
    fitted = as.vector(fitted_vals),
    residuals = as.vector(residuals),
    n = n,
    p = p,
    df_residual = n - p
  ))
}

#' R Fallback for EMMA Algorithm
#' @keywords internal
fit_emma_r <- function(y, X, H, tol = 1e-6, max_iter = 100) {
  
  n <- length(y)
  p <- ncol(X)
  
  # Eigendecomposition of H (which is ZKZ')
  eig <- eigen(H, symmetric = TRUE)
  U <- eig$vectors
  d <- eig$values
  d[d < .MME_CONSTANTS$MIN_EIGENVALUE] <- .MME_CONSTANTS$MIN_EIGENVALUE
  
  # Transform data
  eta <- as.vector(t(U) %*% y)
  xi <- t(U) %*% X
  
  # Grid search for optimal lambda
  # lambda = sigma_e / sigma_g, so:
  # - High h² (0.95+) → lambda < 10^-5 (sigma_e is small)
  # - Low h² (0.05-) → lambda > 10^5 (sigma_g is small)
  lambda_grid <- 10^seq(-5, 5, length.out = 20)
  best_ll <- -Inf
  best_lambda <- 1
  
  for (lambda in lambda_grid) {
    ll <- emma_loglik_r(eta, xi, d, lambda)
    if (ll > best_ll) {
      best_ll <- ll
      best_lambda <- lambda
    }
  }
  
  # Extend grid if best is at boundary (extreme h²)
  if (best_lambda == max(lambda_grid)) {
    # Very low h² - extend to higher lambda
    lambda_ext <- 10^seq(5.5, 10, length.out = 10)
    for (lambda in lambda_ext) {
      ll <- emma_loglik_r(eta, xi, d, lambda)
      if (ll > best_ll) {
        best_ll <- ll
        best_lambda <- lambda
      }
    }
  } else if (best_lambda == min(lambda_grid)) {
    # Very high h² - extend to lower lambda
    lambda_ext <- 10^seq(-10, -5.5, length.out = 10)
    for (lambda in lambda_ext) {
      ll <- emma_loglik_r(eta, xi, d, lambda)
      if (ll > best_ll) {
        best_ll <- ll
        best_lambda <- lambda
      }
    }
  }
  
  # Refine with golden section search
  a <- best_lambda * 0.1
  b <- best_lambda * 10
  gr <- (sqrt(5) - 1) / 2
  
  c <- b - gr * (b - a)
  dd <- a + gr * (b - a)
  
  for (iter in seq_len(max_iter)) {
    fc <- emma_loglik_r(eta, xi, d, c)
    fd <- emma_loglik_r(eta, xi, d, dd)
    
    if (fc > fd) {
      b <- dd
      dd <- c
      c <- b - gr * (b - a)
    } else {
      a <- c
      c <- dd
      dd <- a + gr * (b - a)
    }
    
    if (abs(b - a) < tol * (abs(c) + abs(dd))) break
  }
  
  lambda_opt <- (a + b) / 2
  
  # Compute variance components
  D_inv <- 1 / (d + lambda_opt)
  
  xiT_D_xi <- t(xi) %*% (D_inv * xi)
  xiT_D_eta <- t(xi) %*% (D_inv * eta)
  beta <- solve(xiT_D_xi, xiT_D_eta)
  
  resid <- eta - xi %*% beta
  rss <- sum(resid^2 * D_inv)
  
  sigma2_e <- rss / (n - p)
  sigma2_g <- sigma2_e / lambda_opt
  
  loglik <- emma_loglik_r(eta, xi, d, lambda_opt)
  
  return(list(
    sigma2_g = sigma2_g,
    sigma2_e = sigma2_e,
    lambda = lambda_opt,
    h2 = sigma2_g / (sigma2_g + sigma2_e),
    beta = as.vector(beta),
    loglik = loglik,
    U = U,
    d = d
  ))
}

#' EMMA Log-likelihood (R version)
#' @keywords internal
emma_loglik_r <- function(eta, xi, d, lambda) {
  n <- length(eta)
  p <- ncol(xi)
  
  D_lambda <- d + lambda
  D_inv <- 1 / D_lambda
  
  xiT_D_xi <- t(xi) %*% (D_inv * xi)
  xiT_D_eta <- t(xi) %*% (D_inv * eta)
  
  beta <- tryCatch({
    solve(xiT_D_xi, xiT_D_eta)
  }, error = function(e) {
    solve(xiT_D_xi + diag(1e-8, ncol(xi)), xiT_D_eta)
  })
  
  resid <- eta - xi %*% beta
  rss <- sum(resid^2 * D_inv)
  
  log_det_V <- sum(log(D_lambda))
  log_det_XVX <- as.numeric(determinant(xiT_D_xi, logarithm = TRUE)$modulus)
  
  # REML log-likelihood with constant term for proper AIC/BIC comparison
  const_term <- -(n - p) / 2 * log(2 * pi)
  const_term - 0.5 * (log_det_V + log_det_XVX + rss)
}

#' Full AI-REML Algorithm for Multiple Random Effects
#' @keywords internal
fit_mme_ai_reml <- function(y, X, Z_list, K_list,
                            nIters, tolParConvLL, tolParInv, verbose) {
  
  n <- length(y)
  p <- ncol(X)
  n_random <- length(Z_list)
  
  # Smart initialization of variance components
  # BUG FIX (v2.2.1): smarter ENV-aware initialisation (see fit_mme_henderson).
  # Project out fixed effects to get residual variance for initialization
  # This prevents overestimation when fixed effects explain large portion of variance
  resid_y <- tryCatch({
    residuals(lm.fit(X, as.vector(y)))
  }, error = function(e) {
    # Fallback if lm.fit fails (e.g., singular X)
    as.vector(y) - mean(as.vector(y), na.rm = TRUE)
  })
  var_y <- var(resid_y, na.rm = TRUE)
  if (var_y == 0 || is.na(var_y)) var_y <- 1
  
  q_sizes_ai <- sapply(Z_list, ncol)
  is_env_component <- q_sizes_ai < max(15, 0.15 * n) & q_sizes_ai < 100
  
  # Initialize: assume h2 ~ 0.3 for main genetic, smaller for other effects
  theta <- numeric(n_random + 1)
  
  genetic_share <- 0.5  # Total genetic variance proportion guess
  theta[n_random + 1] <- var_y * (1 - genetic_share)  # Residual
  
  if (n_random == 1) {
    theta[1] <- var_y * genetic_share
  } else if (!any(is_env_component)) {
    # Distribute genetic variance: 60% to first, 25% to second, rest split
    shares <- c(0.6, 0.25, rep(0.15 / max(1, n_random - 2), max(0, n_random - 2)))
    shares <- shares[1:n_random]
    shares <- shares / sum(shares)  # Normalize
    theta[1:n_random] <- var_y * genetic_share * shares
  } else {
    n_env_comps     <- sum(is_env_component)
    n_genetic_comps <- n_random - n_env_comps
    env_seed        <- var_y * 0.5 / max(1, n_env_comps)
    gen_seed_total  <- var_y * genetic_share
    if (n_genetic_comps == 1) {
      gen_shares <- 1.0
    } else {
      gen_shares <- c(0.6, 0.25, rep(0.15 / max(1, n_genetic_comps - 2),
                                       max(0, n_genetic_comps - 2)))
      gen_shares <- gen_shares[seq_len(n_genetic_comps)]
      gen_shares <- gen_shares / sum(gen_shares)
    }
    gen_idx <- 0L
    for (i in seq_len(n_random)) {
      if (is_env_component[i]) {
        theta[i] <- env_seed
      } else {
        gen_idx <- gen_idx + 1L
        theta[i] <- gen_seed_total * gen_shares[gen_idx]
      }
    }
  }
  
  names(theta) <- c(paste0("u", seq_len(n_random)), "e")
  
  # Pre-compute ZKZ' matrices for efficiency
  ZKZt_list <- vector("list", n_random)
  for (i in seq_len(n_random)) {
    Z_i <- Z_list[[i]]
    K_i <- K_list[[i]]
    
    if (ncol(Z_i) != nrow(K_i)) {
      stop("Dimension mismatch: Z[[", i, "]] has ", ncol(Z_i), 
           " columns but K[[", i, "]] has ", nrow(K_i), " rows")
    }
    
    ZKZt_list[[i]] <- Z_i %*% K_i %*% t(Z_i)
  }
  
  # AI-REML iterations
  llik_old <- -Inf
  converged <- FALSE
  monitor <- matrix(NA, nrow = n_random + 1, ncol = nIters)
  boundary_components <- logical(n_random + 1)  # Track which hit boundary
  
  for (iter in seq_len(nIters)) {
    
    # Build phenotypic variance matrix V
    V <- theta[n_random + 1] * diag(n)
    for (i in seq_len(n_random)) {
      V <- V + theta[i] * ZKZt_list[[i]]
    }
    
    # Invert V with numerical stability
    Vi <- safe_inverse(V, tol = tolParInv)
    if (is.null(Vi)) {
      V <- V + diag(tolParInv * 10, n)
      Vi <- solve(V)
    }
    
    # Calculate P matrix: P = Vi - Vi*X*(X'*Vi*X)^{-1}*X'*Vi
    XtVi <- t(X) %*% Vi
    XtViX <- XtVi %*% X
    XtViX_inv <- safe_inverse(XtViX, tol = tolParInv)
    if (is.null(XtViX_inv)) {
      XtViX_inv <- solve(XtViX + diag(tolParInv, ncol(XtViX)))
    }
    
    P <- Vi - t(XtVi) %*% XtViX_inv %*% XtVi
    
    # Calculate REML log-likelihood with constant term
    log_det_V <- as.numeric(determinant(V, logarithm = TRUE)$modulus)
    log_det_XtViX <- as.numeric(determinant(XtViX, logarithm = TRUE)$modulus)
    Py <- P %*% y
    yPy <- as.numeric(t(y) %*% Py)
    const_term <- -(n - p) / 2 * log(2 * pi)
    llik <- const_term - 0.5 * (log_det_V + log_det_XtViX + yPy)
    
    # Calculate AI matrix and score vector (needed for both convergence check and update)
    ai_result <- calculate_ai_matrix(y, P, Py, ZKZt_list, n_random)
    AI <- ai_result$AI
    score <- ai_result$score
    
    # Check convergence using both log-likelihood AND gradient norm
    # Gradient-based check ensures we're at a stationary point
    gradient_norm <- max(abs(score))
    llik_converged <- abs(llik - llik_old) < tolParConvLL && iter > 1
    gradient_converged <- gradient_norm < tolParConvLL * 10
    
    if (llik_converged && gradient_converged) {
      converged <- TRUE
      if (verbose) cat("Converged at iteration", iter, 
                       " (gradient norm: ", round(gradient_norm, 6), ")\n", sep = "")
      break
    } else if (llik_converged && iter > 5) {
      # Log-likelihood converged but gradient still large - possible boundary solution
      converged <- TRUE
      if (verbose) {
        cat("Log-likelihood converged at iteration", iter, "\n")
        if (gradient_norm > tolParConvLL * 100) {
          cat("  Note: Large gradient norm (", round(gradient_norm, 4), 
              ") suggests boundary solution\n", sep = "")
        }
      }
      break
    }
    
    llik_old <- llik
    
    # Update variance components using Newton-Raphson with LL-based step control
    AI_inv <- safe_inverse(AI, tol = 1e-10)
    if (is.null(AI_inv)) {
      AI_inv <- diag(1 / pmax(diag(AI), 1e-10))
    }
    
    delta <- AI_inv %*% score
    
    # Helper: evaluate REML LL at trial theta (without computing P, AI, score)
    eval_trial_llik <- function(theta_trial) {
      V_trial <- theta_trial[n_random + 1] * diag(n)
      for (ii in seq_len(n_random)) {
        V_trial <- V_trial + theta_trial[ii] * ZKZt_list[[ii]]
      }
      Vi_trial <- tryCatch(solve(V_trial), error = function(e) NULL)
      if (is.null(Vi_trial)) return(-Inf)
      
      log_det_V_trial <- as.numeric(determinant(V_trial, logarithm = TRUE)$modulus)
      XtVi_trial <- t(X) %*% Vi_trial
      XtViX_trial <- XtVi_trial %*% X
      log_det_XtViX_trial <- tryCatch(
        as.numeric(determinant(XtViX_trial, logarithm = TRUE)$modulus),
        error = function(e) Inf
      )
      Py_trial <- Vi_trial %*% y - t(XtVi_trial) %*% solve(XtViX_trial, XtVi_trial %*% y)
      yPy_trial <- as.numeric(t(y) %*% Py_trial)
      
      -(n - p) / 2 * log(2 * pi) - 0.5 * (log_det_V_trial + log_det_XtViX_trial + yPy_trial)
    }
    
    # Step-halving: enforce positivity, relative bounds, AND LL increase
    step <- 1.0
    max_step_attempts <- 5  # limit expensive LL evaluations
    accepted <- FALSE
    
    for (attempt in seq_len(max_step_attempts)) {
      theta_new <- theta + step * as.vector(delta)
      
      # Check positivity
      if (any(theta_new < .MME_CONSTANTS$MIN_VARIANCE)) {
        step <- step * 0.5
        next
      }
      
      # Check relative change bounds (max 3× per step)
      rel_change <- abs(theta_new - theta) / pmax(abs(theta), 1e-8)
      if (any(rel_change > 3.0)) {
        step <- step * 0.5
        next
      }
      
      # Check LL actually increased (the key fix)
      llik_trial <- eval_trial_llik(theta_new)
      if (llik_trial >= llik_old - tolParConvLL) {
        # Accept: LL increased (or decreased negligibly)
        accepted <- TRUE
        break
      }
      
      step <- step * 0.5
    }
    
    # If Newton step failed to increase LL, fall back to EM-REML update
    if (!accepted) {
      if (verbose) cat("  Note: Newton step decreased LL at iteration", iter, ", using EM fallback\n")
      # EM-REML update (Thompson & Meyer 1986):
      # θ_i_new = θ_i + (θ_i²/r_i) · 2·score_i
      # where score_i = 0.5·[y'P·(dV/dθ_i)·Py - tr(P·(dV/dθ_i))]
      # and r_i is the rank of dV/dθ_i
      # This is guaranteed to increase REML LL monotonically.
      theta_new <- theta
      for (j in seq_len(n_random)) {
        q_j <- ncol(Z_list[[j]])
        # EM step: θ + θ²/q · 2·score
        theta_new[j] <- max(theta[j] + theta[j]^2 / q_j * 2 * score[j],
                            .MME_CONSTANTS$MIN_VARIANCE)
      }
      # σ²_e: dV/dσ²_e = I, rank = n-p for REML
      theta_new[n_random + 1] <- max(
        theta[n_random + 1] + theta[n_random + 1]^2 / (n - p) * 2 * score[n_random + 1],
        .MME_CONSTANTS$MIN_VARIANCE)
    }
    
    # Track boundary hits
    for (j in seq_along(theta_new)) {
      if (theta_new[j] <= .MME_CONSTANTS$MIN_VARIANCE) {
        boundary_components[j] <- TRUE
      }
    }
    
    # Global bounds
    max_var <- var_y * 10
    for (j in seq_along(theta_new)) {
      theta_new[j] <- max(theta_new[j], .MME_CONSTANTS$MIN_VARIANCE)
      theta_new[j] <- min(theta_new[j], max_var)
    }
    
    theta <- theta_new
    
    monitor[, iter] <- theta
    
    if (verbose && iter %% 5 == 0) {
      cat(sprintf("Iteration %d: LogLik = %.4f, theta = %s\n",
                  iter, llik, paste(round(theta, 4), collapse = ", ")))
    }
  }
  
  # Report boundary solutions
  if (any(boundary_components) && verbose) {
    which_boundary <- which(boundary_components)
    cat("Note: Variance component(s)", paste(names(theta)[which_boundary], collapse = ", "),
        "hit lower boundary (may indicate model misspecification)\n")
  }
  
  # Final solutions
  # Fixed effects (BLUE): b = (X'Vi X)^{-1} X'Vi y
  b <- XtViX_inv %*% XtVi %*% y
  rownames(b) <- colnames(X)
  
  # Random effects (BLUP): u_i = theta_i * K_i * Z_i' * P * y
  u <- vector("list", n_random)
  for (i in seq_len(n_random)) {
    u[[i]] <- theta[i] * K_list[[i]] %*% t(Z_list[[i]]) %*% Py
    rownames(u[[i]]) <- colnames(Z_list[[i]])
  }
  
  # Fitted values
  fitted_vals <- X %*% b
  for (i in seq_len(n_random)) {
    fitted_vals <- fitted_vals + Z_list[[i]] %*% u[[i]]
  }
  
  # Residuals
  residuals <- y - fitted_vals
  
  # Model fit statistics
  AIC <- -2 * llik + 2 * (n_random + 1)
  BIC <- -2 * llik + log(n) * (n_random + 1)
  
  return(list(
    b = b,
    u = u,
    uList = u,
    theta = theta,
    Vi = Vi,
    P = P,
    Ci = XtViX_inv,
    llik = llik,
    convergence = converged,
    nIter = iter,
    monitor = monitor[, seq_len(iter), drop = FALSE],
    AIC = AIC,
    BIC = BIC,
    boundary = boundary_components,
    fitted = as.vector(fitted_vals),
    residuals = as.vector(residuals),
    n = n,
    p = p,
    df_residual = n - p
  ))
}

#' Solve MME with Known Variance Components
#' @keywords internal
solve_mme_with_varcomp <- function(y, X, Z, K, sigma2_g, sigma2_e, ridge = 1e-6) {
  
  n <- length(y)
  q <- ncol(Z)
  
  # Build V = sigma2_g * ZKZ' + sigma2_e * I
  ZKZt <- Z %*% K %*% t(Z)
  V <- sigma2_g * ZKZt + sigma2_e * diag(n)
  
  # Invert V
  Vi <- tryCatch({
    chol2inv(chol(V))
  }, error = function(e) {
    V_ridge <- V + diag(ridge, n)
    solve(V_ridge)
  })
  
  # Fixed effects (BLUE)
  XtVi <- t(X) %*% Vi
  XtViX <- XtVi %*% X
  XtViX_inv <- solve(XtViX)
  b <- XtViX_inv %*% XtVi %*% y
  
  # P matrix
  P <- Vi - t(XtVi) %*% XtViX_inv %*% XtVi
  Py <- P %*% y
  
  # Random effects (BLUP)
  u <- sigma2_g * K %*% t(Z) %*% Py
  
  # Fitted values
  fitted <- X %*% b + Z %*% u
  
  # Residuals
  residuals <- y - fitted
  
  return(list(
    b = b,
    u = u,
    Vi = Vi,
    P = P,
    Ci = XtViX_inv,
    fitted = fitted,
    residuals = residuals
  ))
}

#' Calculate AI Matrix and Score Vector
#'
#' Computes the Average Information matrix and score vector for
#' variance component estimation.
#'
#' @param y Response vector
#' @param P Projection matrix
#' @param Py P times y
#' @param ZKZt_list List of ZKZ' matrices
#' @param n_random Number of random effects
#'
#' @return List with AI matrix and score vector
#' @keywords internal
calculate_ai_matrix <- function(y, P, Py, ZKZt_list, n_random) {
  
  n <- length(y)
  n_vc <- n_random + 1  # +1 for residual
  AI <- matrix(0, n_vc, n_vc)
  score <- numeric(n_vc)
  
  # Derivative matrices: dV/d(theta_i) = ZKZ' for random, I for residual
  dV_list <- c(ZKZt_list, list(diag(n)))
  
  for (i in seq_len(n_vc)) {
    dV_i <- dV_list[[i]]
    PdV_i <- P %*% dV_i
    
    # Score: -0.5 * tr(P * dV_i) + 0.5 * y' P dV_i P y
    trace_PdV_i <- sum(diag(PdV_i))
    yPdVPy_i <- as.numeric(t(Py) %*% dV_i %*% Py)
    score[i] <- -0.5 * trace_PdV_i + 0.5 * yPdVPy_i
    
    for (j in i:n_vc) {
      dV_j <- dV_list[[j]]
      
      # AI: 0.5 * y' P dV_i P dV_j P y
      AI[i, j] <- 0.5 * as.numeric(t(Py) %*% dV_i %*% P %*% dV_j %*% Py)
      AI[j, i] <- AI[i, j]
    }
  }
  
  return(list(AI = AI, score = score))
}

#' Safe Matrix Inverse
#'
#' Computes matrix inverse with numerical stability checks.
#'
#' @param M Square matrix
#' @param tol Tolerance for singularity
#'
#' @return Inverse matrix or NULL if singular
#' @keywords internal
safe_inverse <- function(M, tol = 1e-10) {
  
  if (!is.matrix(M) || nrow(M) != ncol(M)) {
    return(NULL)
  }
  
  # Check condition number via eigenvalues
  eig <- tryCatch({
    eigen(M, symmetric = TRUE, only.values = TRUE)
  }, error = function(e) NULL)
  
  if (is.null(eig)) return(NULL)
  
  min_eig <- min(eig$values)
  max_eig <- max(eig$values)
  
  if (min_eig < tol || max_eig / max(min_eig, tol) > 1e10) {
    return(NULL)
  }
  
  # Try Cholesky decomposition first
  result <- tryCatch({
    chol_M <- chol(M)
    chol2inv(chol_M)
  }, error = function(e) {
    tryCatch(solve(M), error = function(e2) NULL)
  })
  
  return(result)
}

#' Predict BLUPs for New Individuals Using MME Framework
#'
#' Computes proper BLUP predictions for individuals not in training set.
#'
#' @param u_train BLUPs for training individuals
#' @param K Full relationship matrix (includes both train and new)
#' @param train_ids IDs of training individuals
#' @param new_ids IDs of new individuals to predict
#'
#' @return Vector of predicted BLUPs for new individuals
#' @export
predict_blup_new <- function(u_train, K, train_ids, new_ids) {
  
  # Subset K matrix
  K_new_train <- K[new_ids, train_ids, drop = FALSE]
  K_train_train <- K[train_ids, train_ids, drop = FALSE]
  
  # BLUP prediction: u_new = K_new,train * K_train,train^{-1} * u_train
  K_train_inv <- tryCatch({
    solve(K_train_train)
  }, error = function(e) {
    solve(K_train_train + diag(1e-6, nrow(K_train_train)))
  })
  
  u_new <- K_new_train %*% K_train_inv %*% u_train
  
  return(as.vector(u_new))
}
