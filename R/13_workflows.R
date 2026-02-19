#' High-Level Workflow Functions
#'
#' Convenience functions for complete genomic prediction workflows including
#' two-stage analysis, reaction norm models, mixed kernels, cross-validation
#' batteries, and result compilation.
#'
#' @name workflows
#' @keywords internal
NULL


# =============================================================================
# TWO-STAGE BLUES
# =============================================================================

#' Compute Two-Stage BLUEs
#'
#' Performs Stage-1 within-environment analysis to extract tester-adjusted BLUEs
#' (Best Linear Unbiased Estimates), then creates a new GBLUPData object for
#' Stage-2 genomic prediction. This absorbs confounding effects (tester, block,
#' replication) that would otherwise bias GxE variance components.
#'
#' Stage 1 fits \code{trait ~ GID + confounders} within each environment using OLS,
#' then extracts per-GID marginal means. Stage 2 uses these BLUEs as the response
#' in genomic prediction models.
#'
#' @param gdata GBLUPData object with raw observations. Phenotypes must contain
#'   GID, ENV, and trait columns.
#' @param fixed Stage-1 confounders to absorb. A one-sided formula referencing
#'   columns in \code{gdata@@phenotypes} (e.g. \code{~ TESTER}, \code{~ TESTER + BLOCK}).
#'   Default is NULL (no confounders, BLUEs = GID means per ENV).
#' @param verbose Print diagnostics (default TRUE).
#'
#' @return A new GBLUPData object where \code{trait} contains per-GID x ENV BLUEs.
#'   Original genotypes and MAF are preserved.
#' @export
#'
#' @examples
#' \dontrun{
#' gdata_2stage <- compute_two_stage_blues(gdata, fixed = ~ TESTER)
#' model_3 <- fit_model3_AE(gdata_2stage)
#' }
compute_two_stage_blues <- function(gdata, fixed = NULL, verbose = TRUE) {

  if (!inherits(gdata, "GBLUPData")) stop("gdata must be a GBLUPData object")

  pheno <- gdata@phenotypes
  required <- c("GID", "ENV", "trait")
  missing <- setdiff(required, colnames(pheno))
  if (length(missing) > 0) stop("Phenotype data missing columns: ", paste(missing, collapse = ", "))

  # Parse confounders from formula
  confounder_vars <- if (!is.null(fixed)) all.vars(fixed) else character(0)
  for (v in confounder_vars) {
    if (!(v %in% colnames(pheno))) {
      stop("Variable '", v, "' in fixed formula not found in phenotypes. ",
           "Available: ", paste(colnames(pheno), collapse = ", "))
    }
  }

  envs <- sort(unique(pheno$ENV))
  if (verbose) cat("Two-stage BLUE extraction across", length(envs), "environments\n")

  blue_list <- lapply(envs, function(env_name) {
    env_data <- pheno[pheno$ENV == env_name & !is.na(pheno$trait), , drop = FALSE]
    if (length(unique(env_data$GID)) < 2) return(NULL)

    # Build Stage-1 formula: trait ~ GID + confounders
    if (length(confounder_vars) > 0) {
      # Check each confounder has >1 level in this environment
      usable <- sapply(confounder_vars, function(v) length(unique(env_data[[v]])) > 1)
      conf_terms <- confounder_vars[usable]
      if (length(conf_terms) > 0) {
        fml <- as.formula(paste("trait ~ GID +", paste(conf_terms, collapse = " + ")))
      } else {
        fml <- trait ~ GID
      }
    } else {
      fml <- trait ~ GID
    }

    fit <- tryCatch(lm(fml, data = env_data), error = function(e) lm(trait ~ GID, data = env_data))

    # Extract per-GID BLUEs (marginal means from fitted values)
    env_data$fitted <- predict(fit, newdata = env_data)
    blue_by_gid <- tapply(env_data$fitted, env_data$GID, mean, na.rm = TRUE)

    data.frame(
      GID   = names(blue_by_gid),
      ENV   = env_name,
      trait  = as.numeric(blue_by_gid),
      stringsAsFactors = FALSE
    )
  })

  blues <- do.call(rbind, Filter(Negate(is.null), blue_list))
  rownames(blues) <- NULL

  # Filter to genotyped individuals
  blues <- blues[blues$GID %in% rownames(gdata@genotypes), ]

  if (verbose) {
    cat(sprintf("  BLUEs: %d GID x ENV combinations (%d hybrids, %d envs)\n",
                nrow(blues), length(unique(blues$GID)), length(unique(blues$ENV))))
  }

  new("GBLUPData",
      genotypes  = gdata@genotypes,
      phenotypes = blues,
      maf        = gdata@maf,
      metadata   = c(gdata@metadata, list(two_stage = TRUE,
                                           stage1_fixed = deparse(fixed))))
}


# =============================================================================
# REACTION NORM KERNEL
# =============================================================================

#' Build Finlay-Wilkinson Reaction Norm Kernel
#'
#' Constructs a reaction norm kernel from genotype-specific environmental
#' sensitivity (slope). The kernel captures GxE arising from heterogeneous
#' responses to environmental quality: \code{K_rn[i,j] = b_i * b_j * A[i,j]}
#' where \code{b_i} is genotype i's centred FW regression slope.
#'
#' @param gdata GBLUPData object with phenotypes containing GID, ENV, trait.
#' @param A Additive relationship matrix (used as base for kernel scaling).
#' @param min_env Minimum environments per genotype for slope estimation (default: 3).
#'   Genotypes with fewer environments receive a slope of 0 (average sensitivity).
#' @param verbose Print diagnostics (default TRUE).
#'
#' @return A reaction norm kernel matrix (n_gid x n_gid).
#' @export
#'
#' @references Finlay, K.W. and Wilkinson, G.N. (1963) The analysis of
#'   adaptation in a plant-breeding programme. \emph{Australian Journal of
#'   Agricultural Research}, 14, 742--754.
#'
#' @examples
#' \dontrun{
#' K_rn <- build_reaction_norm_kernel(gdata, A)
#' model_rn <- fit_gblup(gdata, K_matrices = list(A = A, RN = K_rn),
#'                        effects = c("A", "ENV", "RN"), fixed = ~ TESTER)
#' }
build_reaction_norm_kernel <- function(gdata, A, min_env = 3, verbose = TRUE) {

  pheno <- gdata@phenotypes

  # Environmental index: mean trait value per environment
  env_means <- tapply(pheno$trait, pheno$ENV, mean, na.rm = TRUE)
  env_score <- as.numeric(scale(env_means))
  names(env_score) <- names(env_means)

  # Per-GID x ENV means
  pheno_rn <- pheno[!is.na(pheno$trait), ]
  pheno_rn$env_score <- env_score[pheno_rn$ENV]
  pheno_rn <- pheno_rn[!is.na(pheno_rn$env_score), ]

  gid_env_means <- aggregate(trait ~ GID + env_score, data = pheno_rn, FUN = mean, na.rm = TRUE)

  # Estimate per-GID slope
  gid_slopes <- do.call(rbind, lapply(split(gid_env_means, gid_env_means$GID), function(dd) {
    if (nrow(dd) < min_env) {
      return(data.frame(GID = dd$GID[1], slope = 0))  # centre = average sensitivity
    }
    fit <- lm(trait ~ env_score, data = dd)
    data.frame(GID = dd$GID[1], slope = coef(fit)["env_score"])
  }))
  rownames(gid_slopes) <- gid_slopes$GID

  # Build slope vector aligned to A
  gids_in_A <- rownames(A)
  slopes_vec <- rep(0, length(gids_in_A))
  names(slopes_vec) <- gids_in_A
  matched <- intersect(gid_slopes$GID, gids_in_A)
  slopes_vec[matched] <- gid_slopes[matched, "slope"]
  slopes_vec <- slopes_vec - mean(slopes_vec)  # centre

  # K_rn = (b ⊗ b') * A  (element-wise)
  K_rn <- (slopes_vec %o% slopes_vec) * A
  dimnames(K_rn) <- dimnames(A)

  if (verbose) {
    cat(sprintf("  Reaction norm kernel: %d x %d, slope range [%.2f, %.2f]\n",
                nrow(K_rn), ncol(K_rn), min(slopes_vec), max(slopes_vec)))
  }

  K_rn
}


# =============================================================================
# MIXED (COMBINED) KERNEL
# =============================================================================

#' Build Trace-Scaled Mixed Kernel
#'
#' Combines two or more relationship matrices using variance-component-weighted,
#' trace-normalised combination. This ensures matrices on different scales
#' contribute proportionally to the combined genetic covariance.
#'
#' The weight for matrix K_i is: \code{w_i = sigma2_i * tr(K_i) / sum(sigma2_j * tr(K_j))}
#'
#' @param K_list Named list of relationship matrices to combine.
#' @param varcomp Named numeric vector of variance components (names matching K_list).
#'   If NULL, equal weights are used.
#' @param ridge Ridge to add to the combined matrix (default: 0.01).
#' @param verbose Print diagnostics (default TRUE).
#'
#' @return A single combined relationship matrix.
#' @export
#'
#' @examples
#' \dontrun{
#' K_comb <- build_mixed_kernel(list(A = A, AA = AA),
#'                               varcomp = c(A = sig2_A, AA = sig2_AA))
#' }
build_mixed_kernel <- function(K_list, varcomp = NULL, ridge = 0.01, verbose = TRUE) {

  if (!is.list(K_list) || length(K_list) < 2) {
    stop("K_list must be a named list with >= 2 matrices")
  }

  n <- nrow(K_list[[1]])
  nms <- names(K_list)

  # Compute traces
  traces <- sapply(K_list, function(K) sum(diag(K)))

  # Normalise each matrix so mean diagonal ≈ 1
  K_norm <- lapply(seq_along(K_list), function(i) K_list[[i]] / (traces[i] / n))

  # Compute weights
  if (is.null(varcomp)) {
    weights <- rep(1 / length(K_list), length(K_list))
    names(weights) <- nms
  } else {
    vc <- varcomp[nms]
    if (any(is.na(vc))) stop("varcomp names must match K_list names")
    contributions <- vc * traces
    weights <- contributions / sum(contributions)
  }

  # Combine
  K_comb <- matrix(0, n, n)
  for (i in seq_along(K_norm)) {
    K_comb <- K_comb + weights[i] * K_norm[[i]]
  }
  dimnames(K_comb) <- dimnames(K_list[[1]])

  # Ridge
  K_comb <- K_comb + diag(ridge, n)

  if (verbose) {
    cat("  Mixed kernel weights:", paste(sprintf("%s=%.3f", nms, weights), collapse = ", "), "\n")
  }

  K_comb
}


# =============================================================================
# MODEL 3R AND 5C PRESETS
# =============================================================================

#' @rdname model_presets
#' @export
fit_model3R_reaction_norm <- function(gdata, A = NULL, ridge = 0.001,
                                       ridge_RN = 0.05, nIters = 500,
                                       verbose = TRUE, fixed = NULL) {
  if (verbose) cat("=== Fitting Model 3R: A + ENV + RN (Reaction Norm) ===\n")

  if (is.null(A)) {
    K <- calc_relationship_matrices(gdata, matrices = "A", ridge = ridge)
    A <- K$A
  }

  K_rn <- build_reaction_norm_kernel(gdata, A, verbose = verbose)

  model <- fit_gblup(gdata,
                     K_matrices = list(A = A, RN = K_rn),
                     effects = c("A", "ENV", "RN"),
                     ridge = ridge,
                     ridge_per_matrix = c(A = ridge, RN = ridge_RN),
                     nIters = nIters, verbose = verbose, fixed = fixed,
                     use.emma = FALSE)
  model@description <- "Model 3R: Reaction Norm (A + ENV + RN)"
  return(model)
}

#' @rdname model_presets
#' @export
fit_model5C_mixed_kernel <- function(gdata, A = NULL, AA = NULL,
                                      varcomp = NULL,
                                      ridge = 0.001, ridge_AE = 0.05,
                                      nIters = 500, verbose = TRUE,
                                      fixed = NULL) {
  if (verbose) cat("=== Fitting Model 5C: K_comb + ENV + AE (Mixed Kernel) ===\n")

  if (is.null(A) || is.null(AA)) {
    K <- calc_relationship_matrices(gdata, matrices = c("A", "AA"), ridge = 0)
    if (is.null(A)) A <- K$A
    if (is.null(AA)) AA <- K$AA
  }

  K_comb <- build_mixed_kernel(list(A = A, AA = AA), varcomp = varcomp,
                                verbose = verbose)

  # Pass A separately (unridged) for AE auto-computation
  model <- fit_gblup(gdata,
                     K_matrices = list(K_comb = K_comb, A = A),
                     effects = c("K_comb", "ENV", "AE"),
                     ridge = ridge,
                     ridge_per_matrix = c(K_comb = 0, AE = ridge_AE),  # K_comb already ridged
                     nIters = nIters, verbose = verbose, fixed = fixed,
                     use.emma = FALSE,
                     K_already_regularized = FALSE)
  model@description <- "Model 5C: Mixed Kernel (K_comb + ENV + AE)"
  return(model)
}


# =============================================================================
# CV BATTERY
# =============================================================================

#' Run Cross-Validation Battery
#'
#' Runs a standardised set of cross-validation schemes across multiple model
#' configurations. Returns a structured data frame of results for direct
#' publication.
#'
#' @param gdata GBLUPData object (raw or 2-stage BLUEs).
#' @param K_matrices List of relationship matrices.
#' @param models Named list of model configurations. Each element is a list with
#'   \code{effects} (character vector) and optionally \code{fixed}, \code{gdata}
#'   (to override per-model), \code{ridge_per_matrix}.
#'   Default: runs M1 through M5 with standard settings.
#' @param schemes Character vector of CV schemes: \code{"CV1"} (random 5-fold),
#'   \code{"CV2"} (leave-one-env-out). Default: \code{c("CV1", "CV2")}.
#' @param n_folds Number of folds for CV1 (default: 5).
#' @param n_reps Number of repetitions (default: 3).
#' @param ridge Base ridge (default: 0.001).
#' @param ridge_per_matrix Named vector of per-matrix ridge values.
#' @param verbose Print progress (default TRUE).
#'
#' @return Data frame with columns: Model, Scheme, Rep, Fold, PA (predictive ability).
#' @export
#'
#' @examples
#' \dontrun{
#' cv_results <- run_cv_battery(gdata, K,
#'                               models = list(M1 = list(effects = "A")))
#' }
run_cv_battery <- function(gdata, K_matrices,
                           models = NULL,
                           schemes = c("CV1", "CV2"),
                           n_folds = 5, n_reps = 3,
                           ridge = 0.001,
                           ridge_per_matrix = NULL,
                           verbose = TRUE) {

  results <- list()

  for (model_name in names(models)) {
    cfg <- models[[model_name]]
    g <- if (!is.null(cfg$gdata)) cfg$gdata else gdata
    eff <- cfg$effects
    fx <- cfg$fixed
    rpm <- if (!is.null(cfg$ridge_per_matrix)) cfg$ridge_per_matrix else ridge_per_matrix

    for (scheme in schemes) {
      if (verbose) cat(sprintf("  CV: %s / %s ...\n", model_name, scheme))

      cv_res <- tryCatch(
        cv_gblup(g, K_matrices, effects = eff,
                 fixed = fx, scheme = scheme,
                 n_folds = n_folds, n_reps = n_reps,
                 ridge = ridge, ridge_per_matrix = rpm,
                 recompute_K = FALSE, verbose = FALSE),
        error = function(e) {
          if (verbose) cat("    FAILED:", conditionMessage(e), "\n")
          NULL
        }
      )

      if (!is.null(cv_res)) {
        pa_vals <- cv_res@metrics$predictive_ability
        results[[length(results) + 1]] <- data.frame(
          Model  = model_name,
          Scheme = scheme,
          PA_mean = mean(pa_vals, na.rm = TRUE),
          PA_sd   = sd(pa_vals, na.rm = TRUE),
          n_folds = length(pa_vals),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  do.call(rbind, results)
}


# =============================================================================
# COMPILE RESULTS
# =============================================================================

#' Compile Results from Multiple Models
#'
#' Extracts variance components, heritability, and fit statistics from a named
#' list of fitted models into a publication-ready summary table.
#'
#' @param models Named list of GBLUPModel objects.
#' @param cv_results Optional data frame from \code{run_cv_battery}.
#'
#' @return A list with elements:
#'   \describe{
#'     \item{varcomp}{Data frame of all variance components across models}
#'     \item{summary}{Data frame with h2, AIC, BIC, convergence per model}
#'     \item{cv}{CV results (if provided)}
#'   }
#' @export
compile_results <- function(models, cv_results = NULL) {

  # Variance components
  vc_list <- lapply(names(models), function(nm) {
    vc <- get_varcomp(models[[nm]])
    vc$Model <- nm
    vc
  })
  vc_all <- do.call(rbind, vc_list)

  # Summary statistics
  summary_df <- data.frame(
    Model = names(models),
    n_VC = sapply(models, function(m) nrow(m@varcomp)),
    Total_Var = sapply(models, function(m) sum(m@varcomp$Variance)),
    Residual = sapply(models, function(m) {
      idx <- grep("Residual", m@varcomp$Component, ignore.case = TRUE)
      if (length(idx) > 0) m@varcomp$Variance[idx] else NA
    }),
    h2 = sapply(models, function(m) {
      h <- tryCatch(heritability(m), error = function(e) c(h2 = NA))
      as.numeric(h["h2"])
    }),
    LogLik = sapply(models, function(m) if (!is.null(m@model$llik)) m@model$llik else NA),
    AIC = sapply(models, function(m) if (!is.null(m@model$AIC)) m@model$AIC else NA),
    BIC = sapply(models, function(m) if (!is.null(m@model$BIC)) m@model$BIC else NA),
    Converged = sapply(models, function(m) {
      if (!is.null(m@model$convergence)) m@model$convergence else NA
    }),
    stringsAsFactors = FALSE
  )

  list(varcomp = vc_all, summary = summary_df, cv = cv_results)
}
