#' Environmental Similarity Kernels for Covariate-Specific G×E
#'
#' Functions to construct environmental similarity matrices (Omega) from
#' user-defined covariates. These matrices replace the identity kernel
#' in standard G×E models, enabling decomposition of G×E variance into
#' covariate-specific components.
#'
#' @name env_kernels
#' @keywords internal
NULL

#' Build Gaussian Environmental Similarity Kernel
#'
#' Constructs an e × e environmental similarity matrix from any numeric
#' covariate vector using a Gaussian (RBF) kernel. Suitable for use with
#' the \code{env_kernels} parameter of \code{fit_gblup()}.
#'
#' @details
#' The Gaussian kernel computes similarity as:
#'
#' \deqn{\Omega_{jj'} = \exp\left(-\frac{(z_j - z_{j'})^2}{2h^2}\right)}{Omega[j,j'] = exp(-(z_j - z_j')^2 / (2h^2))}
#'
#' where \eqn{z_j} is the standardised covariate value for environment j and h is
#' the bandwidth parameter.
#'
#' When \code{bandwidth = NULL} (default), the median heuristic is used:
#' h = median of all pairwise distances. This is a widely-used data-driven
#' bandwidth selection method.
#'
#' @param covariate_values Named numeric vector of covariate values per environment.
#'   Names must be environment IDs matching the ENV column in phenotype data.
#'   The covariate can be any numeric environmental variable: temperature,
#'   rainfall, soil pH, photoperiod, nitrogen application rate, etc.
#' @param bandwidth Optional bandwidth parameter (scalar > 0). Controls the
#'   width of the similarity function: larger h = broader similarity.
#'   If NULL (default), the median heuristic is applied.
#' @return Symmetric positive semi-definite matrix (e × e) with row/column
#'   names matching the input environment names. Diagonal elements are 1.
#'   Off-diagonal elements range from 0 (very different) to 1 (identical).
#' @export
#' @examples
#' \dontrun{
#' # Example 1: Temperature-based environmental kernel
#' temp <- setNames(weather$fill_temp, weather$ENV)
#' Omega_temp <- build_env_kernel(temp)
#'
#' # Example 2: Soil pH kernel
#' soil_ph <- setNames(soil_data$pH, soil_data$location)
#' Omega_soil <- build_env_kernel(soil_ph)
#'
#' # Example 3: Custom bandwidth
#' Omega_rain <- build_env_kernel(rainfall_values, bandwidth = 0.5)
#'
#' # Use in fit_gblup for covariate-specific G×E
#' model <- fit_gblup(
#'   gdata = gdata,
#'   K_matrices = list(A = G),
#'   effects = c("A", "AE_temperature", "AE_soil_ph"),
#'   env_kernels = list(
#'     AE_temperature = Omega_temp,
#'     AE_soil_ph = Omega_soil
#'   )
#' )
#'
#' # Variance decomposition: which covariate shapes G×E?
#' model@varcomp
#' }
build_env_kernel <- function(covariate_values, bandwidth = NULL) {

  # ── Input validation ──
  if (!is.numeric(covariate_values)) {
    stop("covariate_values must be a numeric vector")
  }

  if (is.null(names(covariate_values))) {
    stop("covariate_values must be a named vector (names = environment IDs).\n",
         "  Example: setNames(weather$temperature, weather$ENV)")
  }

  if (length(covariate_values) < 2) {
    stop("At least 2 environments required to build a similarity kernel")
  }

  # ── Handle NAs ──
  n_na <- sum(is.na(covariate_values))
  if (n_na > 0) {
    if (n_na == length(covariate_values)) {
      stop("All covariate values are NA")
    }
    message(sprintf("  %d NA value(s) imputed with covariate mean", n_na))
    covariate_values[is.na(covariate_values)] <- mean(covariate_values, na.rm = TRUE)
  }

  # ── Standardise ──
  z <- (covariate_values - mean(covariate_values)) / sd(covariate_values)

  # ── Distance matrix ──
  D <- as.matrix(dist(z))

  # ── Bandwidth selection ──
  if (is.null(bandwidth)) {
    bandwidth <- median(D[upper.tri(D)], na.rm = TRUE)
    if (is.na(bandwidth) || bandwidth <= 0) bandwidth <- 1
  } else {
    if (!is.numeric(bandwidth) || length(bandwidth) != 1 || bandwidth <= 0) {
      stop("bandwidth must be a positive scalar")
    }
  }

  # ── Gaussian kernel ──
  Omega <- exp(-D^2 / (2 * bandwidth^2))
  rownames(Omega) <- colnames(Omega) <- names(covariate_values)

  # ── Verify positive semi-definiteness ──
  eigenvalues <- eigen(Omega, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigenvalues < -1e-10)) {
    warning("Environmental kernel has negative eigenvalues. ",
            "Adding small ridge to ensure PSD.", call. = FALSE)
    Omega <- Omega + diag(abs(min(eigenvalues)) + 1e-6, nrow(Omega))
  }

  return(Omega)
}

#' Build Linear Environmental Similarity Kernel
#'
#' Constructs an e × e environmental similarity matrix using a linear kernel
#' (dot product of standardised covariates). Suitable when the relationship
#' between environmental similarity and the covariate is expected to be linear.
#'
#' @param covariate_values Named numeric vector of covariate values per environment
#' @return Symmetric positive semi-definite matrix (e × e)
#' @export
#' @examples
#' \dontrun{
#' temp <- setNames(weather$temperature, weather$ENV)
#' Omega_linear <- build_env_kernel_linear(temp)
#' }
build_env_kernel_linear <- function(covariate_values) {

  if (!is.numeric(covariate_values) || is.null(names(covariate_values))) {
    stop("covariate_values must be a named numeric vector")
  }

  # Handle NAs
  if (any(is.na(covariate_values))) {
    covariate_values[is.na(covariate_values)] <- mean(covariate_values, na.rm = TRUE)
  }

  # Standardise
  z <- (covariate_values - mean(covariate_values)) / sd(covariate_values)

  # Linear kernel: Omega = z * z'
  Omega <- outer(z, z)
  rownames(Omega) <- colnames(Omega) <- names(covariate_values)

  return(Omega)
}

#' Build Multi-Covariate Environmental Kernel
#'
#' Constructs an e × e environmental similarity matrix from multiple
#' covariates simultaneously, using a multivariate Gaussian kernel.
#'
#' @param covariate_matrix Matrix or data frame of covariates (rows = environments,
#'   columns = covariates). Row names must be environment IDs.
#' @param bandwidth Optional bandwidth (scalar or per-covariate vector)
#' @return Symmetric positive semi-definite matrix (e × e)
#' @export
#' @examples
#' \dontrun{
#' # Combine multiple covariates into one kernel
#' cov_mat <- data.frame(
#'   temp = weather$fill_temp,
#'   vpd = weather$mean_vpd,
#'   row.names = weather$ENV
#' )
#' Omega_combined <- build_env_kernel_multi(cov_mat)
#' }
build_env_kernel_multi <- function(covariate_matrix, bandwidth = NULL) {

  if (is.data.frame(covariate_matrix)) {
    covariate_matrix <- as.matrix(covariate_matrix)
  }

  if (is.null(rownames(covariate_matrix))) {
    stop("covariate_matrix must have row names (environment IDs)")
  }

  # Standardise each column
  for (j in seq_len(ncol(covariate_matrix))) {
    vals <- covariate_matrix[, j]
    vals[is.na(vals)] <- mean(vals, na.rm = TRUE)
    covariate_matrix[, j] <- (vals - mean(vals)) / sd(vals)
  }

  # Multivariate distance
  D <- as.matrix(dist(covariate_matrix))

  # Bandwidth
  if (is.null(bandwidth)) {
    bandwidth <- median(D[upper.tri(D)], na.rm = TRUE)
    if (is.na(bandwidth) || bandwidth <= 0) bandwidth <- 1
  }

  Omega <- exp(-D^2 / (2 * bandwidth^2))
  rownames(Omega) <- colnames(Omega) <- rownames(covariate_matrix)

  return(Omega)
}
