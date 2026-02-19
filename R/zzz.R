#' @useDynLib CompGEBLUP, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

.onAttach <- function(libname, pkgname) {
  pkg_version <- tryCatch(
    as.character(utils::packageVersion(pkgname)),
    error = function(e) "1.0"
  )
  
  packageStartupMessage(
    "CompGEBLUP v", pkg_version, "\n",
    "Comprehensive Genomic Prediction with GBLUP\n",
    "Author: Muhammad Kamran\n",
    "Type vignette('user_guide', package='CompGEBLUP') for help\n\n",
    "Core Models (preset functions):\n",
    "  Model 1: A                     fit_model1_A()\n",
    "  Model 2: A + D                 fit_model2_AD()\n",
    "  Model 3: A + ENV + AE          fit_model3_AE()\n",
    "  Model 4: A + D + ENV + AE      fit_model4_ADE()      # DE optional\n",
    "  Model 5: A + ENV + AE + AA     fit_model5_ADAA()      # D optional\n\n",
    "Custom effects supported: pass any named kernel in K_matrices\n",
    "  e.g. effects = c('A', 'ENV', 'RN') with K_matrices = list(A=..., RN=...)\n"
  )
}

# Package-level options
.CompGEBLUP_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  # Check for OpenMP support (set once at load time)
  .CompGEBLUP_env$openmp_available <- check_openmp_support()
  invisible(NULL)
}

#' Check OpenMP Support
#' 
#' @return Logical indicating if OpenMP is available
#' @keywords internal
check_openmp_support <- function() {
  tryCatch({
    TRUE
  }, error = function(e) FALSE)
}

#' Get Package Options
#' 
#' @return List of package options
#' @keywords internal
get_compgeblup_options <- function() {
  list(
    openmp_available = .CompGEBLUP_env$openmp_available
  )
}
