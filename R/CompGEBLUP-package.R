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
#'   \item Henderson EM-REML with full trace-corrected variance component updates
#'   \item AI-REML with Aitken acceleration for multi-random effect models
#'   \item Exact REML log-likelihood for reliable convergence monitoring
#'   \item Support for additive, dominance, and epistatic effects
#'   \item Genotype-by-environment interaction modeling
#'   \item Custom effect support: arbitrary user-defined relationship matrices
#'   \item Per-matrix differential ridge regularisation (\code{ridge_per_matrix})
#'   \item Multi-trait genomic selection
#'   \item Cross-validation without data leakage (CV0, CV1, CV2 schemes)
#'   \item EMMAX-style fast GWAS with pre-computed decomposition
#'   \item Fast C++ matrix operations via RcppArmadillo
#'   \item Optional OpenMP parallelization for GWAS
#' }
#'
#' @section Core Models:
#' The 5 core models have preset functions for quick fitting. Model 4 drops DE
#' by default (overfitting risk), and Model 5 drops D by default (boundary
#' estimate when epistasis is present). Both can be re-enabled via toggles.
#' \itemize{
#'   \item Model 1: A (Baseline additive) — \code{\link{fit_model1_A}}
#'   \item Model 2: A + D (Dominance) — \code{\link{fit_model2_AD}}
#'   \item Model 3: A + ENV + AE (Additive GxE) — \code{\link{fit_model3_AE}}
#'   \item Model 4: A + D + ENV + AE (GxE, DE optional) — \code{\link{fit_model4_ADE}}
#'   \item Model 5: A + ENV + AE + AA (Epistasis, D optional) — \code{\link{fit_model5_ADAA}}
#' }
#'
#' @section Custom Effects:
#' Beyond the 7 built-in effect types (A, D, ENV, AE, DE, AA, AAE), users can
#' pass any named relationship matrix in \code{K_matrices} and reference it by
#' name in \code{effects}. The package auto-detects whether the kernel is
#' GID-level or GID.ENV-level based on dimensions. This enables reaction norm
#' kernels, combined kernels, RKHS kernels, or any user-defined covariance
#' structure. Example:
#' \preformatted{
#' fit_gblup(gdata,
#'           K_matrices = list(A = A, RN = K_reaction_norm),
#'           effects = c("A", "ENV", "RN"))
#' }
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{simulate_genotypes}}: Generate simulated genotype data
#'   \item \code{\link{simulate_phenotypes}}: Generate simulated phenotypes
#'   \item \code{\link{calc_relationship_matrices}}: Calculate relationship matrices
#'   \item \code{\link{fit_gblup}}: Fit GBLUP model (supports \code{ridge_per_matrix})
#'   \item \code{\link{cv_gblup}}: Cross-validation
#'   \item \code{\link{gwas}}: Genome-wide association study
#'   \item \code{\link{fit_multi_trait}}: Multi-trait model
#'   \item \code{\link{make_positive_definite}}: Ridge or eigenvalue-floor regularisation
#'   \item \code{\link{add_ridge}}: Diagonal ridge regularisation
#' }
#'
#' @section C++ Functions:
#' \itemize{
#'   \item \code{\link{calc_A_matrix_cpp}}: Fast additive relationship matrix
#'   \item \code{\link{calc_D_matrix_cpp}}: Fast dominance relationship matrix
#'   \item \code{\link{emma_reml_cpp}}: Efficient REML variance components
#'   \item \code{\link{gwas_emmax_cpp}}: Fast GWAS marker testing
#'   \item \code{\link{fit_henderson_mme_cpp}}: C++ Henderson MME solver
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
