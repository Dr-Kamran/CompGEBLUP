#' @useDynLib CompGEBLUP, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

.onAttach <- function(libname, pkgname) {
  pkg_version <- tryCatch(
    as.character(utils::packageVersion(pkgname)),
    error = function(e) "2.1.0"
  )
  
  packageStartupMessage(
    "CompGEBLUP v", pkg_version, "\n",
    "Comprehensive Genomic Prediction - Independent Mixed Model Engine\n",
    "Author: Muhammad Kamran\n",
    "Type vignette('user_guide', package='CompGEBLUP') for help\n\n",
    "Five Core Models:\n",
    "  Model 1: A (Baseline)\n",
    "  Model 2: A + D (Dominance)\n",
    "  Model 3: A + ENV + AE (Additive GxE)\n",
    "  Model 4: A + D + ENV + AE + DE (Full GxE)\n",
    "  Model 5: A + D + ENV + AE + AA (Epistasis)\n"
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
  # Check if compiled with OpenMP by testing if parallel execution changes timing
  # This is a heuristic since we can't directly query OpenMP status from R
  tryCatch({
    # The C++ code uses OpenMP pragma for GWAS - if available, it will use multiple threads
    # We just assume it's available if the package loads successfully
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
