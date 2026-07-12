#' Calculate Relationship Matrices for GBLUP
#'
#' Wrapper function to calculate multiple relationship matrices at once
#'
#' @param gdata GBLUPData object containing genotypes
#' @param matrices Character vector specifying which matrices to calculate.
#'   Options: "A" (additive), "D" (dominance), "AA" (additive epistasis),
#'   "AE" (additive x environment), "DE" (dominance x environment),
#'   "AAE" (epistatic x environment), "AD" (additive x dominance), "DD" (dominance x dominance)
#' @param method Method for calculating A matrix ("VanRaden", "UAR", "UARadj")
#' @param min.MAF Minimum minor allele frequency for filtering
#' @param ridge Ridge parameter for numerical stability
#'
#' @return Named list of relationship matrices
#' @export
#' @examples
#' \dontrun{
#' gdata <- simulate_genotypes(n_ind = 100, n_snp = 500)
#' gdata <- simulate_phenotypes(gdata, n_env = 3, h2 = 0.6)
#' 
#' # Calculate just A matrix
#' K <- calc_relationship_matrices(gdata, matrices = "A")
#' 
#' # Calculate multiple matrices
#' K <- calc_relationship_matrices(gdata, matrices = c("A", "D", "AE"))
#' }
calc_relationship_matrices <- function(gdata, 
                                       matrices = "A",
                                       method = c("VanRaden", "UAR", "UARadj"),
                                       min.MAF = 0.01,
                                       ridge = 0.001) {
  
  # Match arguments
  method <- match.arg(method)
  
  if (!inherits(gdata, "GBLUPData")) {
    stop("gdata must be a GBLUPData object")
  }
  
  # Validate matrices argument
  valid_matrices <- c("A", "D", "AA", "AD", "DD", "AE", "DE", "AAE")
  invalid <- setdiff(matrices, valid_matrices)
  if (length(invalid) > 0) {
    stop("Invalid matrix types: ", paste(invalid, collapse = ", "),
         ". Valid options: ", paste(valid_matrices, collapse = ", "))
  }
  
  # Get genotype matrix
  M <- gdata@genotypes
  
  if (is.null(M) || nrow(M) == 0) {
    stop("No genotype data available in gdata")
  }
  
  # Get number of environments if needed for GxE matrices
  n_env <- length(unique(gdata@phenotypes$ENV))
  env_names <- unique(as.character(gdata@phenotypes$ENV))
  
  # Initialize result list
  K_list <- list()
  
  # Calculate each requested matrix
  for (mat_type in matrices) {
    
    if (mat_type == "A") {
      K_list$A <- build_A_matrix(M, method = method, min.MAF = min.MAF)
      
    } else if (mat_type == "D") {
      K_list$D <- build_D_matrix(M, min.MAF = min.MAF)
      
    } else if (mat_type == "AA") {
      if (!"A" %in% names(K_list)) {
        A <- build_A_matrix(M, method = method, min.MAF = min.MAF)
      } else {
        A <- K_list$A
      }
      K_list$AA <- build_E_matrix(M, type = "A#A", A = A, min.MAF = min.MAF)
      
    } else if (mat_type == "AD") {
      if (!"A" %in% names(K_list)) {
        A <- build_A_matrix(M, method = method, min.MAF = min.MAF)
      } else {
        A <- K_list$A
      }
      if (!"D" %in% names(K_list)) {
        D <- build_D_matrix(M, min.MAF = min.MAF)
      } else {
        D <- K_list$D
      }
      K_list$AD <- build_E_matrix(M, type = "A#D", A = A, D = D, min.MAF = min.MAF)
      
    } else if (mat_type == "DD") {
      if (!"D" %in% names(K_list)) {
        D <- build_D_matrix(M, min.MAF = min.MAF)
      } else {
        D <- K_list$D
      }
      K_list$DD <- build_E_matrix(M, type = "D#D", D = D, min.MAF = min.MAF)
      
    } else if (mat_type == "AE") {
      if (!"A" %in% names(K_list)) {
        A <- build_A_matrix(M, method = method, min.MAF = min.MAF)
      } else {
        A <- K_list$A
      }
      K_list$AE <- build_GE_matrix(A, n_env = n_env, env_names = env_names)
      
    } else if (mat_type == "DE") {
      if (!"D" %in% names(K_list)) {
        D <- build_D_matrix(M, min.MAF = min.MAF)
      } else {
        D <- K_list$D
      }
      K_list$DE <- build_GE_matrix(D, n_env = n_env, env_names = env_names)
      
    } else if (mat_type == "AAE") {
      if (!"AA" %in% names(K_list)) {
        if (!"A" %in% names(K_list)) {
          A <- build_A_matrix(M, method = method, min.MAF = min.MAF)
        } else {
          A <- K_list$A
        }
        AA <- build_E_matrix(M, type = "A#A", A = A, min.MAF = min.MAF)
      } else {
        AA <- K_list$AA
      }
      K_list$AAE <- build_GE_matrix(AA, n_env = n_env, env_names = env_names)
    }
  }
  
  # Add ridge to all matrices for numerical stability
  if (ridge > 0) {
    K_list <- lapply(K_list, function(K) {
      K + diag(ridge, nrow(K), ncol(K))
    })
  }
  
  return(K_list)
}

#' Extract Genotype Matrix from GBLUPData
#'
#' @param gdata GBLUPData object
#' @param individuals Optional vector of individual IDs to extract
#'
#' @return Genotype matrix
#' @export
#' @examples
#' \dontrun{
#' gdata <- simulate_genotypes(n_ind = 100, n_snp = 500)
#' M <- get_genotypes(gdata)
#' M_subset <- get_genotypes(gdata, individuals = c("ID1", "ID2", "ID3"))
#' }
get_genotypes <- function(gdata, individuals = NULL) {
  
  if (!inherits(gdata, "GBLUPData")) {
    stop("gdata must be a GBLUPData object")
  }
  
  M <- gdata@genotypes
  
  if (!is.null(individuals)) {
    # Extract specific individuals
    missing <- individuals[!individuals %in% rownames(M)]
    if (length(missing) > 0) {
      stop("Individuals not found in genotype matrix: ", 
           paste(head(missing, 5), collapse = ", "),
           if (length(missing) > 5) paste0(" (and ", length(missing) - 5, " more)"))
    }
    M <- M[individuals, , drop = FALSE]
  }
  
  return(M)
}

#' Extract Phenotype Data from GBLUPData
#'
#' @param gdata GBLUPData object
#' @param trait Optional trait name to extract
#' @param environment Optional environment to filter
#'
#' @return Data frame of phenotypes
#' @export
#' @examples
#' \dontrun{
#' gdata <- simulate_genotypes(n_ind = 100, n_snp = 500)
#' gdata <- simulate_phenotypes(gdata, n_env = 3, h2 = 0.6)
#' pheno <- get_phenotypes(gdata)
#' pheno_e1 <- get_phenotypes(gdata, environment = "E1")
#' }
get_phenotypes <- function(gdata, trait = NULL, environment = NULL) {
  
  if (!inherits(gdata, "GBLUPData")) {
    stop("gdata must be a GBLUPData object")
  }
  
  pheno <- gdata@phenotypes
  
  if (!is.null(environment)) {
    if (!environment %in% pheno$ENV) {
      stop("Environment not found: ", environment,
           ". Available: ", paste(unique(pheno$ENV), collapse = ", "))
    }
    pheno <- pheno[pheno$ENV == environment, ]
  }
  
  if (!is.null(trait)) {
    if (!trait %in% names(pheno)) {
      stop("Trait not found in phenotype data: ", trait,
           ". Available: ", paste(setdiff(names(pheno), c("GID", "ENV")), collapse = ", "))
    }
    # Return only relevant columns
    pheno <- pheno[, c("GID", "ENV", trait)]
  }
  
  return(pheno)
}

#' Validate Relationship Matrix
#'
#' Checks if a matrix is a valid relationship matrix
#'
#' @param K Matrix to validate
#' @param check_pd Check if positive definite (default TRUE)
#' @param tol Tolerance for checks
#'
#' @return TRUE if valid, otherwise throws error
#' @export
#' @examples
#' \dontrun{
#' K <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
#' validate_K_matrix(K)
#' }
validate_K_matrix <- function(K, check_pd = TRUE, tol = 1e-6) {
  
  # Check if matrix
  if (!is.matrix(K)) {
    stop("K must be a matrix")
  }
  
  # Check if square
  if (nrow(K) != ncol(K)) {
    stop("K must be a square matrix")
  }
  
  # Check if symmetric
  if (!isSymmetric(K, tol = tol)) {
    warning("K is not symmetric, symmetrizing...")
    K <- (K + t(K)) / 2
  }
  
  # Check if has dimnames
  if (is.null(rownames(K)) || is.null(colnames(K))) {
    warning("K should have dimnames (individual IDs)")
  }
  
  # Check if positive definite
  if (check_pd) {
    eig <- eigen(K, symmetric = TRUE, only.values = TRUE)
    min_eig <- min(eig$values)
    
    if (min_eig < -tol) {
      stop("K is not positive semi-definite. Minimum eigenvalue: ", round(min_eig, 6))
    }
    
    if (min_eig < tol) {
      warning("K is near singular. Minimum eigenvalue: ", round(min_eig, 6),
              ". Consider adding a ridge parameter.")
    }
  }
  
  return(TRUE)
}

#' Align Matrices to Common Individuals
#'
#' Ensures all matrices and data have the same individuals in the same order
#'
#' @param K_list List of relationship matrices
#' @param phenotypes Phenotype data frame
#' @param id_col Name of ID column in phenotypes (default "GID")
#'
#' @return List with aligned matrices and phenotypes
#' @export
#' @examples
#' \dontrun{
#' K <- calc_relationship_matrices(gdata, matrices = "A")
#' aligned <- align_matrices(K, gdata@phenotypes)
#' }
align_matrices <- function(K_list, phenotypes, id_col = "GID") {
  
  if (!is.list(K_list)) {
    stop("K_list must be a list of matrices")
  }
  
  if (!id_col %in% names(phenotypes)) {
    stop("id_col '", id_col, "' not found in phenotypes")
  }
  
  # Get IDs from phenotypes
  pheno_ids <- unique(as.character(phenotypes[[id_col]]))
  
  # Get IDs from each matrix
  matrix_ids <- lapply(K_list, rownames)
  
  # Find common IDs
  common_ids <- pheno_ids
  for (ids in matrix_ids) {
    if (!is.null(ids)) {
      common_ids <- intersect(common_ids, ids)
    }
  }
  
  if (length(common_ids) == 0) {
    stop("No common individuals found between phenotypes and relationship matrices")
  }
  
  n_removed <- length(pheno_ids) - length(common_ids)
  if (n_removed > 0) {
    message("Removed ", n_removed, " individuals not present in all matrices")
  }
  
  # Subset and reorder matrices
  K_aligned <- lapply(K_list, function(K) {
    if (!is.null(rownames(K))) {
      K[common_ids, common_ids, drop = FALSE]
    } else {
      K
    }
  })
  
  # Subset phenotypes
  pheno_aligned <- phenotypes[phenotypes[[id_col]] %in% common_ids, ]
  
  return(list(
    K_matrices = K_aligned,
    phenotypes = pheno_aligned,
    common_ids = common_ids
  ))
}

#' Calculate Genomic Prediction Accuracy
#'
#' Computes correlation between observed and predicted values
#'
#' @param observed Vector of observed values
#' @param predicted Vector of predicted values
#' @param method Correlation method ("pearson", "spearman")
#'
#' @return Correlation coefficient
#' @export
#' @examples
#' \dontrun{
#' observed <- c(1, 2, 3, 4, 5)
#' predicted <- c(1.1, 2.2, 2.9, 4.1, 4.8)
#' calc_accuracy(observed, predicted)
#' }
calc_accuracy <- function(observed, predicted, method = c("pearson", "spearman")) {
  
  method <- match.arg(method)
  
  # Remove missing values
  complete <- complete.cases(observed, predicted)
  
  if (sum(complete) < 2) {
    warning("Less than 2 complete observations")
    return(NA)
  }
  
  obs <- observed[complete]
  pred <- predicted[complete]
  
  # Calculate correlation
  acc <- cor(obs, pred, method = method)
  
  return(acc)
}

#' Calculate Mean Squared Error
#'
#' @param observed Vector of observed values
#' @param predicted Vector of predicted values
#'
#' @return MSE value
#' @export
#' @examples
#' \dontrun{
#' observed <- c(1, 2, 3, 4, 5)
#' predicted <- c(1.1, 2.2, 2.9, 4.1, 4.8)
#' calc_mse(observed, predicted)
#' }
calc_mse <- function(observed, predicted) {
  
  # Remove missing values
  complete <- complete.cases(observed, predicted)
  
  if (sum(complete) == 0) {
    warning("No complete observations")
    return(NA)
  }
  
  obs <- observed[complete]
  pred <- predicted[complete]
  
  mse <- mean((obs - pred)^2)
  
  return(mse)
}

#' Calculate Root Mean Squared Error
#'
#' @param observed Vector of observed values
#' @param predicted Vector of predicted values
#'
#' @return RMSE value
#' @export
calc_rmse <- function(observed, predicted) {
  sqrt(calc_mse(observed, predicted))
}

#' Calculate Mean Absolute Error
#'
#' @param observed Vector of observed values
#' @param predicted Vector of predicted values
#'
#' @return MAE value
#' @export
calc_mae <- function(observed, predicted) {
  
  complete <- complete.cases(observed, predicted)
  
  if (sum(complete) == 0) {
    warning("No complete observations")
    return(NA)
  }
  
  mean(abs(observed[complete] - predicted[complete]))
}