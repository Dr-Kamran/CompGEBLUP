#' CompGEBLUP: Comprehensive Genomic Prediction with Multi-Trait GBLUP
#'
#' @description
#' A complete toolkit for genomic selection including additive, dominance,
#' epistatic effects, GxE interactions, multi-trait selection, cross-validation,
#' and GWAS. Features an independent REML-based mixed model engine with
#' EMMA-style efficient algorithms and fast C++ matrix operations using RcppArmadillo.
#'
#' @section Key Features:
#' \itemize{
#'   \item Independent mixed model engine (no external dependencies)
#'   \item EMMA-style efficient REML for O(n) per-iteration complexity
#'   \item AI-REML for multi-random effect models
#'   \item Support for additive, dominance, and epistatic effects
#'   \item Genotype-by-environment interaction modeling
#'   \item Multi-trait genomic selection
#'   \item Cross-validation without data leakage (CV0, CV1, CV2 schemes)
#'   \item EMMAX-style fast GWAS with pre-computed decomposition
#'   \item Fast C++ matrix operations via RcppArmadillo
#'   \item Optional OpenMP parallelization for GWAS
#' }
#'
#' @section Core Models:
#' \itemize{
#'   \item Model 1: A (Baseline additive)
#'   \item Model 2: A + D (Dominance)
#'   \item Model 3: A + ENV + AE (Additive GxE)
#'   \item Model 4: A + D + ENV + AE + DE (Full GxE)
#'   \item Model 5: A + D + ENV + AE + AA (Epistasis)
#' }
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{simulate_genotypes}}: Generate simulated genotype data
#'   \item \code{\link{simulate_phenotypes}}: Generate simulated phenotypes
#'   \item \code{\link{calc_relationship_matrices}}: Calculate relationship matrices
#'   \item \code{\link{fit_gblup}}: Fit GBLUP model
#'   \item \code{\link{cv_gblup}}: Cross-validation
#'   \item \code{\link{gwas}}: Genome-wide association study
#'   \item \code{\link{fit_multi_trait}}: Multi-trait model
#' }
#'
#' @section C++ Functions:
#' \itemize{
#'   \item \code{\link{calc_A_matrix_cpp}}: Fast additive relationship matrix
#'   \item \code{\link{calc_D_matrix_cpp}}: Fast dominance relationship matrix
#'   \item \code{\link{emma_reml_cpp}}: Efficient REML variance components
#'   \item \code{\link{gwas_emmax_cpp}}: Fast GWAS marker testing
#' }
#'
#' @docType package
#' @name CompGEBLUP-package
#' @aliases CompGEBLUP
#' @useDynLib CompGEBLUP, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Matrix Diagonal sparseMatrix crossprod tcrossprod solve t
#' @importFrom methods new setClass setMethod setGeneric show
#' @importFrom stats aggregate as.formula coef complete.cases cor cov2cor
#'   density lm median pnorm ppoints qbeta qchisq reshape runif sd var
#'   rnorm rbinom prcomp pt na.pass na.omit pchisq model.matrix
#'   as.dist hclust setNames cov
#' @importFrom grDevices colorRampPalette rainbow rgb hcl.colors adjustcolor
#' @importFrom graphics abline arrows axis barplot boxplot hist image
#'   legend lines par pie polygon segments text plot.new title points rect
#' @importFrom utils head txtProgressBar setTxtProgressBar
#' @keywords internal
"_PACKAGE"
