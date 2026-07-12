#' Build a Location-Recurrent, Climate-Orthogonalized Environmental Kernel (CW-RN)
#'
#' Constructs the structured environmental covariance Omega* used by the
#' structured major-gene G x E model (sm-GEBLUP). It fuses two pieces:
#'   (1) I_L  -- a LOCATION-identity block: environments at the same location
#'       covary, capturing the RECURRENT part of G x E (a site property that
#'       repeats across years);
#'   (2) Omega_clim -- a climate-similarity kernel built ONLY from
#'       between-location climate normals, ORTHOGONALIZED against I_L so it does
#'       not double-count location identity. This carries the part that
#'       transfers to a genuinely NEW location.
#'
#' The weather-deviation part (within-location, year-specific) is deliberately
#' excluded from prediction: it is transient and unforecastable at deployment.
#'
#' @param env_names character vector of environment labels (year-location or
#'   location codes), length n_env, matching the phenotype ENV coding.
#' @param location factor/character of length n_env giving each environment's
#'   LOCATION (so multiple year-locations can share a location).
#' @param clim_features optional numeric matrix (n_env x p) of environmental
#'   covariates; the between-location means become the climate normals. If NULL,
#'   Omega* reduces to the pure recurrent kernel I_L (+ ridge).
#' @param w_clim weight (0..1) on the orthogonalized climate block relative to
#'   the recurrent block. Default 0.5.
#' @param ridge diagonal ridge for positive-definiteness (default 0.05).
#'
#' @return list with Omega (n_env x n_env, the kernel to pass to
#'   build_env_qtl_matrices), plus the diagnostic align_before / align_after
#'   (correlation of the climate kernel with location identity, pre/post
#'   orthogonalization) so leakage can be reported.
#' @export
build_cwrn_kernel <- function(env_names, location, clim_features = NULL,
                              w_clim = 0.5, ridge = 0.05) {
  n <- length(env_names)
  stopifnot(length(location) == n)
  loc <- as.character(location)

  # recurrent block: 1 if same location
  B_L <- outer(loc, loc, "==") * 1.0

  offd <- function(A) A[upper.tri(A)]
  align_before <- NA_real_; align_after <- NA_real_
  Om_clim <- matrix(0, n, n)

  if (!is.null(clim_features)) {
    X <- as.matrix(clim_features)
    for (j in seq_len(ncol(X))) { v <- X[, j]; v[!is.finite(v)] <- mean(v[is.finite(v)], na.rm = TRUE); X[, j] <- v }
    # between-location climate normals (each env gets its location's mean)
    clim <- t(vapply(seq_len(n), function(i) colMeans(X[loc == loc[i], , drop = FALSE]), numeric(ncol(X))))
    Zc <- scale(clim); Zc[!is.finite(Zc)] <- 0
    Om_clim_raw <- tcrossprod(Zc) / ncol(Zc)
    align_before <- suppressWarnings(cor(offd(Om_clim_raw), offd(B_L)))

    # ORTHOGONALIZE the climate kernel against the location-identity block:
    # regress out B_L from Om_clim (both taken as vectors over off-diagonals),
    # so the climate block cannot re-encode location identity.
    b <- offd(B_L); k <- offd(Om_clim_raw)
    beta <- sum((b - mean(b)) * (k - mean(k))) / sum((b - mean(b))^2)
    resid <- k - beta * (b - mean(b))
    Om_clim <- matrix(0, n, n); Om_clim[upper.tri(Om_clim)] <- resid
    Om_clim <- Om_clim + t(Om_clim); diag(Om_clim) <- 1
    align_after <- suppressWarnings(cor(offd(Om_clim), offd(B_L)))
  }

  Omega <- (1 - w_clim) * B_L + (if (!is.null(clim_features)) w_clim * Om_clim else 0)
  diag(Omega) <- diag(Omega) + ridge
  dimnames(Omega) <- list(env_names, env_names)

  ev <- min(eigen(Omega, symmetric = TRUE, only.values = TRUE)$values)
  if (ev <= 1e-8) { Omega <- Omega + diag(1e-6 - ev, n); diag(Omega) <- diag(Omega) }

  list(Omega = Omega, align_before = align_before, align_after = align_after,
       recurrent_block = B_L)
}
