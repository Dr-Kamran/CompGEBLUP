################################################################################
# G2F CompGEBLUP Analysis Script
# 
# Corrected analysis workflow for Genomes to Fields (G2F) data using the
# CompGEBLUP package. This script demonstrates proper genomic prediction
# with additive, dominance, and epistatic relationship matrices.
#
# Key Improvements Over Standalone Scripts:
# 1. Uses package's fast C++ PLINK BED reader (decode_bed_cpp)
# 2. Uses package's AI-REML for robust variance component estimation
# 3. Proper hybrid D matrix calculation using Vitezica et al. (2013)
# 4. Cross-validation for prediction accuracy assessment
# 5. Phenotype outlier filtering
# 6. Proper multi-environment handling
# 7. Set random seeds for reproducibility
#
# References:
# - VanRaden (2008): Additive relationship matrix
# - Vitezica et al. (2013): Dominance relationship matrix
# - Su et al. (2012): Epistatic relationship matrices
# - Kang et al. (2008): EMMA algorithm for efficient REML
################################################################################

library(CompGEBLUP)

# Set random seed for reproducibility
set.seed(42)

################################################################################
# SECTION 1: Load and Prepare Data
################################################################################

cat("=" , rep("=", 78), "\n", sep = "")
cat("G2F COMPGEBLUP ANALYSIS\n")
cat("=" , rep("=", 78), "\n\n", sep = "")

# ----- 1.1: Load Genotype Data -----
cat("Step 1: Loading genotype data...\n")

# Replace these paths with your actual data files
plink_prefix <- "path/to/g2f_data"  # e.g., "/data/G2F/genotypes"
bed_file <- paste0(plink_prefix, ".bed")
bim_file <- paste0(plink_prefix, ".bim")
fam_file <- paste0(plink_prefix, ".fam")

# Check if files exist
if (!file.exists(bed_file)) {
  stop("BED file not found. Please update the 'plink_prefix' path.")
}

# Load PLINK binary files using fast C++ decoder
# Returns genotypes coded as 0/1/2 (number of alternate alleles)
plink_data <- load_plink_bed(
  bed_file = bed_file,
  bim_file = bim_file,
  fam_file = fam_file,
  coding = "additive",  # Convert to -1/0/1 for CompGEBLUP
  impute_missing = TRUE,
  verbose = TRUE
)

# Extract genotype matrix (individuals x SNPs, coded as -1/0/1)
M_inbreds <- plink_data$geno
cat(sprintf("  Loaded %d inbred lines with %d SNPs\n", 
            nrow(M_inbreds), ncol(M_inbreds)))

# Clean up large object
rm(plink_data)
gc()


# ----- 1.2: Load Phenotype Data -----
cat("\nStep 2: Loading phenotype data...\n")

# Replace with your actual phenotype file
pheno_file <- "path/to/phenotypes.csv"
if (!file.exists(pheno_file)) {
  stop("Phenotype file not found. Please update the 'pheno_file' path.")
}

pheno_raw <- read.csv(pheno_file, stringsAsFactors = FALSE)

# Expected columns: hybrid_id, environment, trait (e.g., Yield_Mg_ha)
# Also need parent IDs for hybrid analysis: parent1, parent2

required_cols <- c("hybrid_id", "environment", "trait", "parent1", "parent2")
if (!all(required_cols %in% colnames(pheno_raw))) {
  stop("Phenotype data must contain columns: ", paste(required_cols, collapse = ", "))
}

cat(sprintf("  Loaded %d phenotype records from %d environments\n",
            nrow(pheno_raw), length(unique(pheno_raw$environment))))


# ----- 1.3: Phenotype Outlier Filtering -----
cat("\nStep 3: Filtering phenotype outliers...\n")

# Remove extreme outliers within each environment (±4 SD)
pheno_filtered <- do.call(rbind, lapply(split(pheno_raw, pheno_raw$environment), function(env_data) {
  trait_mean <- mean(env_data$trait, na.rm = TRUE)
  trait_sd <- sd(env_data$trait, na.rm = TRUE)
  
  # Keep observations within ±4 SD
  lower_bound <- trait_mean - 4 * trait_sd
  upper_bound <- trait_mean + 4 * trait_sd
  
  outliers <- env_data$trait < lower_bound | env_data$trait > upper_bound
  n_outliers <- sum(outliers, na.rm = TRUE)
  
  if (n_outliers > 0) {
    cat(sprintf("  Environment %s: Removed %d outliers (%.1f%%)\n",
                env_data$environment[1], n_outliers, 
                100 * n_outliers / nrow(env_data)))
  }
  
  env_data[!outliers, ]
}))

cat(sprintf("  Retained %d records after outlier filtering\n", nrow(pheno_filtered)))


# ----- 1.4: Construct Hybrid Genotype Matrix -----
cat("\nStep 4: Constructing hybrid genotype matrix...\n")

# Get unique hybrids
hybrids <- unique(pheno_filtered$hybrid_id)
n_hybrids <- length(hybrids)

# Initialize hybrid genotype matrix
M_hybrids <- matrix(NA, nrow = n_hybrids, ncol = ncol(M_inbreds))
rownames(M_hybrids) <- hybrids
colnames(M_hybrids) <- colnames(M_inbreds)

# For each hybrid, compute genotype as mean of parents
# This assumes F1 hybrids: M_hybrid = (M_parent1 + M_parent2) / 2
for (i in 1:n_hybrids) {
  hybrid <- hybrids[i]
  
  # Get parent IDs for this hybrid
  idx <- which(pheno_filtered$hybrid_id == hybrid)[1]
  p1 <- pheno_filtered$parent1[idx]
  p2 <- pheno_filtered$parent2[idx]
  
  # Check if both parents are in genotype data
  if (p1 %in% rownames(M_inbreds) && p2 %in% rownames(M_inbreds)) {
    M_hybrids[i, ] <- (M_inbreds[p1, ] + M_inbreds[p2, ]) / 2
  } else {
    warning(sprintf("Parents %s and/or %s not found for hybrid %s", p1, p2, hybrid))
  }
}

# Remove hybrids with missing parent genotypes
missing_idx <- apply(M_hybrids, 1, function(x) all(is.na(x)))
if (any(missing_idx)) {
  cat(sprintf("  Warning: Removing %d hybrids with missing parent genotypes\n",
              sum(missing_idx)))
  M_hybrids <- M_hybrids[!missing_idx, ]
  hybrids <- hybrids[!missing_idx]
}

cat(sprintf("  Constructed genotypes for %d hybrids\n", nrow(M_hybrids)))

# Clean up
rm(M_inbreds)
gc()


################################################################################
# SECTION 2: Build Relationship Matrices
################################################################################

cat("\n", rep("=", 80), "\n", sep = "")
cat("SECTION 2: Building Relationship Matrices\n")
cat(rep("=", 80), "\n\n", sep = "")

# ----- 2.1: Create GBLUPData Object -----
cat("Step 5: Creating GBLUPData object...\n")

# Prepare phenotype data frame in CompGEBLUP format
pheno_df <- data.frame(
  GID = pheno_filtered$hybrid_id,
  ENV = pheno_filtered$environment,
  trait = pheno_filtered$trait,
  stringsAsFactors = FALSE
)

# Keep only hybrids with genotype data
pheno_df <- pheno_df[pheno_df$GID %in% rownames(M_hybrids), ]

cat(sprintf("  Final dataset: %d observations on %d hybrids in %d environments\n",
            nrow(pheno_df), 
            length(unique(pheno_df$GID)),
            length(unique(pheno_df$ENV))))

# Create GBLUPData object
gdata <- new("GBLUPData",
             genotypes = M_hybrids,
             phenotypes = pheno_df)


# ----- 2.2: Build Additive Relationship Matrix -----
cat("\nStep 6: Building additive relationship matrix (A)...\n")

# Use VanRaden (2008) method with C++ acceleration
A <- build_A_matrix(
  M = M_hybrids,
  method = "VanRaden",
  min.MAF = 0.01,
  use.cpp = TRUE,
  verbose = TRUE
)

cat(sprintf("  A matrix: %d x %d, mean diagonal = %.3f\n",
            nrow(A), ncol(A), mean(diag(A))))


# ----- 2.3: Build Dominance Relationship Matrix -----
cat("\nStep 7: Building dominance relationship matrix (D)...\n")

# Use Vitezica et al. (2013) method with C++ acceleration
# This is the CORRECT method for dominance (not Su et al. 2012)
D <- build_D_matrix(
  M = M_hybrids,
  method = "Vitezica",
  min.MAF = 0.01,
  use.cpp = TRUE,
  verbose = TRUE
)

cat(sprintf("  D matrix: %d x %d, mean diagonal = %.3f\n",
            nrow(D), ncol(D), mean(diag(D))))

# NOTE: The package's build_D_matrix() using Vitezica et al. (2013) correctly
# computes the dominance matrix diagonal from heterozygosity patterns.
# This is the CORRECT method (not Su et al. 2012 used in some standalone scripts).
#
# For reference, standalone scripts sometimes incorrectly used:
#   d_ii <- 1 - A[p1i, p2i]  # WRONG
# The correct approach uses the D matrix itself:
#   d_ii <- max(D[p1i, p2i], 0.01)  # Could be used for manual adjustment
# However, the package function handles this automatically, so no manual
# adjustment is needed.


# ----- 2.4: Build Epistatic Relationship Matrix -----
cat("\nStep 8: Building epistatic relationship matrix (AA)...\n")

# Additive x Additive epistasis (A#A)
AA <- build_E_matrix(
  M = M_hybrids,
  type = "A#A",
  A = A,
  min.MAF = 0.01,
  use.cpp = TRUE
)

cat(sprintf("  AA matrix: %d x %d, mean diagonal = %.3f\n",
            nrow(AA), ncol(AA), mean(diag(AA))))

gc()


################################################################################
# SECTION 3: Variance Component Estimation with AI-REML
################################################################################

cat("\n", rep("=", 80), "\n", sep = "")
cat("SECTION 3: Estimating Variance Components\n")
cat(rep("=", 80), "\n\n", sep = "")

# The package's fit_gblup() uses:
# - EMMA for single random effect (very fast)
# - AI-REML for multiple random effects (robust convergence)
# Both are MUCH better than naive EM with ridge inside the loop

# ----- 3.1: Fit Additive Model (A) -----
cat("Model 1: Additive (A)\n")
K_A <- list(A = A)
model_A <- fit_gblup(
  gdata = gdata,
  K_matrices = K_A,
  effects = "A",
  ridge = 0.001,
  nIters = 500,
  tolParConvLL = 1e-6,
  verbose = TRUE,
  use.emma = TRUE  # Use EMMA for speed with single random effect
)

# Extract variance components
vc_A <- get_varcomp(model_A)
h2_A <- heritability(model_A, n_env = length(unique(pheno_df$ENV)))

cat("\nVariance Components (A model):\n")
print(vc_A)
cat(sprintf("  Heritability: %.3f\n\n", h2_A))


# ----- 3.2: Fit Additive + Dominance Model (A+D) -----
cat("Model 2: Additive + Dominance (A+D)\n")
K_AD <- list(A = A, D = D)
model_AD <- fit_gblup(
  gdata = gdata,
  K_matrices = K_AD,
  effects = c("A", "D"),
  ridge = 0.001,
  nIters = 500,
  tolParConvLL = 1e-6,
  verbose = TRUE,
  use.emma = FALSE  # Use AI-REML for multiple variance components
)

vc_AD <- get_varcomp(model_AD)
h2_AD <- heritability(model_AD, n_env = length(unique(pheno_df$ENV)))

cat("\nVariance Components (A+D model):\n")
print(vc_AD)
cat(sprintf("  Narrow-sense heritability: %.3f\n", h2_AD["h2"]))
cat(sprintf("  Broad-sense heritability: %.3f\n\n", h2_AD["H2"]))


# ----- 3.3: Fit Additive + Dominance + Epistasis Model (A+D+AA) -----
cat("Model 3: Additive + Dominance + Epistasis (A+D+AA)\n")
K_ADAA <- list(A = A, D = D, AA = AA)
model_ADAA <- fit_gblup(
  gdata = gdata,
  K_matrices = K_ADAA,
  effects = c("A", "D", "AA"),
  ridge = 0.001,
  nIters = 500,
  tolParConvLL = 1e-6,
  verbose = TRUE,
  use.emma = FALSE
)

vc_ADAA <- get_varcomp(model_ADAA)
h2_ADAA <- heritability(model_ADAA, n_env = length(unique(pheno_df$ENV)))

cat("\nVariance Components (A+D+AA model):\n")
print(vc_ADAA)
cat(sprintf("  Heritability: %.3f\n\n", h2_ADAA))

gc()


################################################################################
# SECTION 4: Cross-Validation for Prediction Accuracy
################################################################################

cat("\n", rep("=", 80), "\n", sep = "")
cat("SECTION 4: Cross-Validation\n")
cat(rep("=", 80), "\n\n", sep = "")

# Set seed for reproducible fold assignment
set.seed(123)

# ----- 4.1: CV for Model A -----
cat("Cross-validating Model A (Additive)...\n")
cv_A <- cv_gblup(
  gdata = gdata,
  K_matrices = K_A,
  effects = "A",
  scheme = "CV1",  # Random k-fold
  n_folds = 5,
  n_reps = 3,
  ridge = 0.001,
  recompute_K = TRUE,  # Avoid data leakage (slower but correct)
  verbose = TRUE
)

cat("\nModel A Cross-Validation Results:\n")
cat(sprintf("  Mean Predictive Ability: %.3f ± %.3f\n",
            mean(cv_A@metrics$predictive_ability, na.rm = TRUE),
            sd(cv_A@metrics$predictive_ability, na.rm = TRUE)))
cat(sprintf("  Mean MSE: %.3f ± %.3f\n\n",
            mean(cv_A@metrics$mse, na.rm = TRUE),
            sd(cv_A@metrics$mse, na.rm = TRUE)))


# ----- 4.2: CV for Model A+D -----
cat("Cross-validating Model A+D...\n")
cv_AD <- cv_gblup(
  gdata = gdata,
  K_matrices = K_AD,
  effects = c("A", "D"),
  scheme = "CV1",
  n_folds = 5,
  n_reps = 3,
  ridge = 0.001,
  recompute_K = TRUE,
  verbose = TRUE
)

cat("\nModel A+D Cross-Validation Results:\n")
cat(sprintf("  Mean Predictive Ability: %.3f ± %.3f\n",
            mean(cv_AD@metrics$predictive_ability, na.rm = TRUE),
            sd(cv_AD@metrics$predictive_ability, na.rm = TRUE)))
cat(sprintf("  Mean MSE: %.3f ± %.3f\n\n",
            mean(cv_AD@metrics$mse, na.rm = TRUE),
            sd(cv_AD@metrics$mse, na.rm = TRUE)))


# ----- 4.3: CV for Model A+D+AA -----
cat("Cross-validating Model A+D+AA...\n")
cv_ADAA <- cv_gblup(
  gdata = gdata,
  K_matrices = K_ADAA,
  effects = c("A", "D", "AA"),
  scheme = "CV1",
  n_folds = 5,
  n_reps = 3,
  ridge = 0.001,
  recompute_K = TRUE,
  verbose = TRUE
)

cat("\nModel A+D+AA Cross-Validation Results:\n")
cat(sprintf("  Mean Predictive Ability: %.3f ± %.3f\n",
            mean(cv_ADAA@metrics$predictive_ability, na.rm = TRUE),
            sd(cv_ADAA@metrics$predictive_ability, na.rm = TRUE)))
cat(sprintf("  Mean MSE: %.3f ± %.3f\n\n",
            mean(cv_ADAA@metrics$mse, na.rm = TRUE),
            sd(cv_ADAA@metrics$mse, na.rm = TRUE)))


################################################################################
# SECTION 5: Model Comparison and Summary
################################################################################

cat("\n", rep("=", 80), "\n", sep = "")
cat("SECTION 5: Model Comparison\n")
cat(rep("=", 80), "\n\n", sep = "")

# Summary table
model_comparison <- data.frame(
  Model = c("A", "A+D", "A+D+AA"),
  Heritability = c(
    h2_A["h2"],  # Use consistent indexing
    h2_AD["h2"],  # narrow-sense h2
    h2_ADAA["h2"]
  ),
  Predictive_Ability = c(
    mean(cv_A@metrics$predictive_ability, na.rm = TRUE),
    mean(cv_AD@metrics$predictive_ability, na.rm = TRUE),
    mean(cv_ADAA@metrics$predictive_ability, na.rm = TRUE)
  ),
  MSE = c(
    mean(cv_A@metrics$mse, na.rm = TRUE),
    mean(cv_AD@metrics$mse, na.rm = TRUE),
    mean(cv_ADAA@metrics$mse, na.rm = TRUE)
  )
)

cat("Model Comparison Summary:\n")
print(model_comparison, row.names = FALSE)

cat("\n\nAnalysis complete!\n")
cat("=" , rep("=", 78), "\n", sep = "")


################################################################################
# SECTION 6: Optional - Predict GEBVs for New Hybrids
################################################################################

# If you have new hybrids (not in training data), you can predict their GEBVs:
# 
# # Load new hybrid genotypes
# M_new <- ... # matrix of new hybrid genotypes
# 
# # Create newdata GBLUPData object
# newdata <- new("GBLUPData",
#                genotypes = M_new,
#                phenotypes = data.frame(GID = rownames(M_new),
#                                       ENV = "Env1",
#                                       trait = NA))
# 
# # Compute relationship matrix between new and training hybrids
# K_new <- list(A = ... )  # compute using build_A_matrix with combined data
# 
# # Predict
# predictions <- predict_gblup(model_AD, newdata = newdata, K_new = K_new)


################################################################################
# NOTES AND KEY IMPROVEMENTS
################################################################################

# 1. PLINK BED PARSING:
#    - Old: Pure R nested loops (hours for 1000 ind x 500k SNPs)
#    - New: C++ decode_bed_cpp() (seconds)
#
# 2. DOMINANCE MATRIX:
#    - Old: Su et al. (2012) coding with buggy hybrid diagonal d_ii = 1 - A[p1,p2]
#    - New: Vitezica et al. (2013) via build_D_matrix() with correct diagonal
#
# 3. VARIANCE COMPONENT ESTIMATION:
#    - Old: Naive EM with ridge inside loop, max 100 iter, poor convergence
#    - New: EMMA (1 random effect) or AI-REML (multiple effects) via fit_mme()
#           with 500+ iterations, log-likelihood tracking, step-halving
#
# 4. CROSS-VALIDATION:
#    - Old: None
#    - New: k-fold CV with recompute_K to avoid data leakage
#
# 5. OUTLIER FILTERING:
#    - Old: None
#    - New: ±4 SD within each environment
#
# 6. MULTI-ENVIRONMENT:
#    - Old: aggregate() averaging loses G×E information
#    - New: Keep individual records, fit environment as fixed effect
#
# 7. REPRODUCIBILITY:
#    - Old: No set.seed()
#    - New: set.seed() for fold assignment and other stochastic operations
#
# 8. MEMORY MANAGEMENT:
#    - Old: No gc() calls
#    - New: gc() after large matrix operations

################################################################################
