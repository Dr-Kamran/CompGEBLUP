# CompGEBLUP 2.3.1

## Critical Bug Fixes

### Henderson EM-REML σ²_e update was simplified (CRITICAL)
Both `fit_mme_henderson` and `fit_mme_henderson_ai` used `σ²_e = RSS/(n-p)` — the
naive OLS residual estimate, not the full EM-REML update. The correct formula is
`σ²_e = (RSS + σ²_e · tr(C⁻¹ · C_data)) / (n - p)`. The trace correction accounts
for variance "borrowed" from the residual by random effects. Without it, σ²_e is
overestimated and genetic VCs underestimated for multi-effect models. The same fix
was applied to the C++ `fit_henderson_mme_cpp`.

### Approximate REML LL in Henderson EM-REML
`fit_mme_henderson` used `-0.5 * [n·log(σ²_e) + RSS/σ²_e]` — the Gaussian LL for
i.i.d. errors — instead of the exact REML LL. This caused unreliable convergence
monitoring and premature convergence. Now computes exact REML LL with `log|C|`,
`log|G|`, and `y'Py` terms in both R and C++ solvers.

### make_positive_definite `ridge=` silently converted to eigenvalue floor
When called with `ridge = 0.05`, the function converted it to `tol = 0.05` and
floored eigenvalues — a fundamentally different operation that distorts matrix
structure. Now `ridge=` adds diagonal ridge (K + ridge*I) as expected.

### cv_fold used eigenvalue floor instead of diagonal ridge
`cv_fold` with `recompute_K=TRUE` called `make_positive_definite(K, tol=ridge)`
instead of `make_positive_definite(K, ridge=ridge)`. This meant CV folds used
eigenvalue flooring while the full-data model used diagonal ridge.

### AA matrix inconsistency between full-data and CV folds
CV folds computed `AA_train = K*K` (raw Hadamard) while full-data used
`build_E_matrix(type="A#A")`. Now CV folds use `build_E_matrix` with Hadamard fallback.

### tolParInv coupled to user's ridge parameter
`fit_gblup` passed the user's `ridge` value as `tolParInv` to the solver. When
`ridge = 0.05`, this set the solver's internal stability parameter too high.
Now uses a fixed `1e-6`.

### LL plateau convergence too loose
Both Henderson solvers declared convergence when `llik_diff < 0.01` after 15
iterations. Tightened to `0.001` after 25 iterations.

## New Features

### Custom effect support in `build_design_matrices`
Any effect name not in the 7 built-in types (A, D, ENV, AE, DE, AA, AAE) is now
treated as a custom effect. A matching K matrix must exist in `K_matrices`. The
function auto-detects whether the kernel is GID-level (n_gid × n_gid) or
GID.ENV-level (n_gid*n_env × n_gid*n_env) based on dimensions. Unknown effects
without matching matrices now produce a warning instead of being silently dropped.

This enables:
- Reaction norm kernels: `effects = c("A", "ENV", "RN")`
- Combined kernels: `effects = c("K_comb", "ENV", "AE")`
- Factor-analytic kernels, RKHS kernels, or any user-defined structure

### `ridge_per_matrix` in `fit_gblup()`
Differential ridge values can now be specified per matrix type, matching the
existing `cv_gblup()` interface.

### Updated model presets
- `fit_model4_ADE()`: DE dropped by default (`include_DE = FALSE`). DE causes
  severe overfitting in hybrid breeding data.
- `fit_model5_ADAA()`: D dropped by default (`include_D = FALSE`). σ²_D often
  hits boundary when epistasis is in the model.
- All presets use `ridge_per_matrix` for differential regularisation.
- Default `nIters` increased from 50 to 500 (50 was insufficient for multi-effect).
- `fit_model1_A()` uses `use.emma = TRUE`; all multi-effect presets use `FALSE`.

### CV fold custom effect pass-through
`cv_fold` with `recompute_K=TRUE` now passes custom-effect K matrices through
from the original `K_matrices` when no built-in rebuild path exists.

### AE Kronecker documentation
`build_design_matrices` documents that auto-computed AE inherits ridge from A.

## Testing
Added `test-v231-fixes.R` with 11 new tests covering all fixes.

---

# CompGEBLUP 2.2.0

## Bug Fixes

### Bug 1: C++ Convention Mixing in Henderson Solver (CRITICAL)
**File:** `src/matrix_ops.cpp`, function `solve_henderson_mme_cpp()`

Used `lambda * K_inv` where `lambda = sigma2_e/sigma2_g`, while data blocks were
divided by `sigma2_e`. Mixed Convention 1/2 makes G⁻¹ penalty σ²_e times too strong.

**Fixed:** `K_inv / sigma2_g` (consistent Convention 2).
**Affects:** `use.emma = TRUE` (default) with single random effect.

### Bug 2: C++ EM Trace Correction (CRITICAL)  
**File:** `src/matrix_ops.cpp`, function `fit_henderson_mme_cpp()`

EM trace correction used `sigma2_e * trace(C_inv_uu * K_inv)`. In Convention 2,
C⁻¹_uu is PEV directly — no σ²_e scaling needed.

**Fixed:** `trace(C_inv_uu * K_inv)` (no σ²_e multiplier).
**Affects:** Same as Bug 1.

### Bug 3: R AI-REML Log-Likelihood (CRITICAL)
**File:** `R/09_mme_solver.R`, function `fit_mme_henderson_ai()`

Log-likelihood used `rss/σ²_e` instead of Henderson's `y'Py = y'y/σ²_e - rhs'·sol`.
The missing `û'K⁻¹û/σ²_g` terms caused LL to decrease during optimization.

**Fixed:** Proper y'Py via Henderson identity.
**Affects:** `use.emma = FALSE` Henderson AI-REML solver.

### Bug 4: R AI-REML Score for σ²_e (CRITICAL)
**File:** `R/09_mme_solver.R`, function `fit_mme_henderson_ai()`

Two errors in the residual variance score:
(a) Used `∂C/∂σ²_e = -C_data/σ²_e` instead of `-C_data/σ²_e²` — derivative wrong
    by factor σ²_e, making the trace term ~550x too large.
(b) Used `(n-p)/σ²_e` instead of `n/σ²_e` — the fixed-effect adjustment is already
    captured by `∂log|C|/∂σ²_e` in Henderson space.

**Fixed:** Score = `0.5*(-n/σ²_e + (rss + tr(C⁻¹·C_data))/σ²_e²)`.
Fisher Information derivative also corrected to `-C_data/σ²_e²`.
**Affects:** All Henderson AI-REML models.

# CompGEBLUP 2.1.2

## Critical Bug Fixes

### Fixed REML Variance Initialization for Models with Large Fixed Effects (CRITICAL)
**Previous bug:** When users specified `fixed = ~ TESTER` (or any formula with many levels), the EM-REML and AI-REML solvers failed to estimate genetic variance correctly. The variance initialization used `var(y)` (total phenotypic variance), but when fixed effects absorbed most of the phenotypic variation, the initial `genetic_share = 0.5` of total variance was far too high. The optimizer then rapidly drove σ²_g → 0, producing h² = 0.

**Fixed:** Now initialize variance components using the residual variance after projecting out fixed effects via `residuals(lm.fit(X, y))`. This ensures that the starting values are appropriate for the actual genetic and residual variance components, leading to proper convergence even with many fixed effect levels.

**Impact:** Users can now reliably use fixed effects (e.g., tester effects in hybrid trials) without encountering variance component estimation failures.

**Locations fixed:** `fit_mme_ai_reml()`, `fit_mme_henderson()` (both R and C++ paths) in `R/09_mme_solver.R`.

### Bug 1: Fixed `compute_pev()` to Return Per-Individual Values (CRITICAL)
**Previous bug:** `compute_pev()` returned a single scalar value (variance of residuals) repeated for all individuals, making all individuals appear equally reliable.

**Fixed:** Now returns per-individual PEV values based on prediction quality (r² within environment groups). Individuals with higher prediction accuracy get lower PEV, correctly reflecting their reliability.

**Impact:** `compute_reliability()` now also returns meaningful per-individual values instead of identical reliability for everyone.

### Bug 2: Fixed `predict_gblup()` K Matrix Construction (CRITICAL)  
**Previous bug:** Used incorrect approximation `K_new %*% t(K_new) / ncol(K_new)` for K_new,new block, which produced wrong scale and meaning.

**Fixed:** Removed the incorrect approximation. Now properly constructs combined K matrix using only the necessary blocks: K[train,train] and K[new,train], which is all that BLUP prediction requires.

## Other Fixes

### Issue 3: Inline h² Calculation Clarified
Added note in verbose output that the simplified inline h² calculation may differ from `heritability()` function for multi-environment models, directing users to use `heritability()` for proper estimates.

### Issue 6: Parameter Naming Corrected
Renamed `fit_emma_r()` parameter from misleading `K` to `H` to accurately reflect that it receives `ZKZt` (H matrix), not the original K matrix.

### Issue 7: REML Log-Likelihood Constant Term Added
Added constant term `-(n-p)/2 * log(2π)` to REML log-likelihood calculations in both `emma_loglik_r()` and `fit_mme_ai_reml()`. This ensures AIC/BIC values are comparable with other software (ASReml, lme4).

### Issue 8: G2F Example Script Fixed
Fixed `inst/examples/G2F_analysis.R` to use lowercase `mse` instead of `MSE` to match actual column names in CV metrics dataframe.

### Issue 9: Model Presets Now Support Fixed Effects
All model preset functions (`fit_model1_A()` through `fit_model5_ADAA()`) now accept a `fixed` parameter to specify fixed effects (e.g., `fixed = ~ TESTER`). This allows users of the preset functions to include tester or other fixed effects in their models.

## Notes

- Issues 4 & 5 (fixed effects pass-through) were already implemented correctly in the codebase
- All changes maintain backward compatibility except for the improved accuracy of PEV/reliability calculations
- Updated tests to reflect new behavior (per-individual variation, no approximation warnings)

---

# CompGEBLUP 2.1.1

## New Features

### Model Comparison Function
New `compare_models_cv()` function to easily compare multiple effect combinations:
```r
effects_to_compare <- list(
  "Additive" = "A",
  "Add+Dom" = c("A", "D"),
  "Add+Dom+Epi" = c("A", "D", "AA")
)
comparison <- compare_models_cv(gdata, K, effects_to_compare)
```
Returns ranked table with statistical tests for significant differences.

### Improved Heritability Output
`heritability()` now returns BOTH interpretations for multi-environment models:
```r
h <- heritability(model)
h["h2"]       # For selection on mean (G×E as noise)
h["h2_total"] # For variance partitioning (G×E as genetic)
```

### CV Progress Indicators with ETA
Cross-validation now shows estimated time remaining:
```
Fold 3 of 10 - ETA: 2.5 min
```

## Bug Fixes (v2.1.1)

### Error Handling: Proper Re-throwing in GWAS
EMMA/GWAS now properly throws errors if both C++ and R fallbacks fail.

### All-NA Columns Handled in `build_AD_matrices()`
Columns with all NA values are now handled gracefully.

### Environment Correlation Validation in `build_GE_matrix()`
Warns if `env_names` doesn't match `rownames(env_cor)`.

### Integer Overflow Protection
Memory calculations use `as.numeric()` to prevent overflow for large datasets.

### CV0 Stricter Validation
Now requires minimum 10 individuals per environment (was 5).

---

# CompGEBLUP 2.1.0

## Critical Bug Fixes

### Heritability Calculation Formula (FIXED - Major)
**Previous bug:** G×E interaction variance was incorrectly included in the NUMERATOR.
This caused inflated heritability estimates in multi-environment models.

**Correct formula:**
```
h² = Va / (Va + Vae/n_env + Ve)
H² = (Va + Vd) / (Va + Vd + Vae/n_env + Vde/n_env + Ve)
```

**Key insight:** Interaction variance (Vae) is NOISE for selection on mean performance,
so it appears in the denominator but NOT the numerator.

### FDR Implementation (FIXED - Critical)
The Benjamini-Hochberg FDR procedure was incorrectly implemented.

**Previous (wrong):**
```r
fdr_threshold <- max(p_sorted[p_sorted <= k / m * alpha], 0)
```

**Correct:**
```r
significant <- p_sorted <= (k / m) * alpha
max_k <- max(which(significant))
fdr_threshold <- p_sorted[max_k]
```

### EMMA Extreme Heritability (FIXED)
Grid search for lambda now extends for extreme h² cases:
- Very high h² (> 0.95): extends grid to 10^-10
- Very low h² (< 0.05): extends grid to 10^10

### CV0 Unbalanced Data Validation (NEW)
CV0 (within-environment) now validates that each environment has at least 
5 individuals. With fewer, results are unreliable.

### Cross-Validation Data Leakage Warning (IMPROVED)
When using `recompute_K = FALSE`, a warning is now issued:
```
Warning: recompute_K = FALSE: Using pre-computed K matrices may cause 
data leakage (test genotypes influence K). Prediction accuracies may be 
inflated by 5-15%. Set recompute_K = TRUE for unbiased CV (slower but correct).
```

### Improved REML Initialization (FIXED)
The AI-REML algorithm now uses smarter variance component initialization:
- First genetic effect (usually additive) gets ~30% of phenotypic variance
- Additional effects get progressively smaller initial values
- Residual initialized at ~50%
- This improves convergence for models with dominance/epistasis

### Boundary Solution Detection (NEW)
REML now tracks and reports when variance components hit boundaries:
- Warns when components are constrained to minimum value
- Suggests possible model misspecification

### Double Ridge Regularization (FIXED)
Previously, ridge was added twice when using `calc_relationship_matrices()` followed by `fit_gblup()`:
1. In `make_positive_definite()` via eigenvalue flooring
2. In `fit_gblup()` via diagonal addition

**Fix:** Added `K_already_regularized` parameter to `fit_gblup()`:
```r
# If K matrices were already regularized (e.g., via make_positive_definite)
model <- fit_gblup(gdata, K, effects = "A", K_already_regularized = TRUE)
```

Cross-validation now automatically sets `K_already_regularized = TRUE` when it builds K matrices internally.

### Accuracy vs Predictive Ability (CLARIFIED)

**Predictive ability** is now clearly documented as THE primary metric:
- `predictive_ability = cor(GEBV, phenotype)` - this is what most GP papers report
- A value of ~0.4 is typical for small training populations with moderate h²

We removed the automatic "true accuracy" calculation because:
- It required estimated h² from training data, which has high variance
- The formula `accuracy = r / sqrt(h²)` assumes TRUE h², not an estimate
- Using estimated h² produces misleading/inflated values

**For users who know their TRUE h²** (e.g., from simulation):
```r
# Convert predictive ability to true accuracy (only with known true h²)
true_acc <- calc_prediction_accuracy(0.41, h2 = 0.6)  # Returns: 0.53

# Calculate theoretical maximum accuracy (Daetwyler formula)
# Note: Me (effective markers) strongly affects this - typical range 100-1000
expected_accuracy(n_train = 80, h2 = 0.6, Me = 100)  # Returns: 0.57
expected_accuracy(n_train = 80, h2 = 0.6, Me = 500)  # Returns: 0.30
```

**For simulated data**, you can validate against known true breeding values:
```r
# Direct calculation when TBV is available
true_accuracy <- cor(GEBV, true_breeding_value)
```

### True Breeding Values Now Stored in Simulations (NEW)
`simulate_phenotypes()` now stores true breeding values in the output:
- `TBV`: True additive breeding value (for calculating true prediction accuracy)
- `TGV`: True total genetic value (includes dominance, epistasis effects)
- Simulation parameters stored as attribute for reference

New function `calc_true_accuracy()` for simulation studies:
```r
# Simulate data with known TBV
gdata <- simulate_genotypes(n_ind = 200, n_snp = 1000)
gdata <- simulate_phenotypes(gdata, h2 = 0.6)

# Fit model and calculate true accuracy
K <- calc_relationship_matrices(gdata, matrices = "A")
model <- fit_gblup(gdata, K, effects = "A")
acc <- calc_true_accuracy(model, gdata)
# Returns: true_accuracy, predictive_ability, theoretical_ratio
```

## Performance Improvements

### Combined A+D Matrix Computation (NEW)
New function `build_AD_matrices()` computes both additive and dominance matrices
in a single pass, avoiding redundant MAF filtering and imputation:
```r
matrices <- build_AD_matrices(M)
A <- matrices$A
D <- matrices$D
```
For n=5000 individuals with both matrices, this is ~40% faster than calling
`build_A_matrix()` and `build_D_matrix()` separately.

### REML Algorithm Constants Documented
Magic numbers in the MME solver are now fully documented with:
- Scientific rationale
- Literature references (ASReml, WOMBAT, BLUPF90)
- User-adjustment guidance

## New Feature: Publication-Quality Visualizations

Updated and added visualization functions for genomic analysis with publication-quality output:

### Updated/New Visualization Functions

- `plot_manhattan()`: Now supports multi-trait Manhattan plots
  - Pass a list of GWAS results for multi-trait visualization
  - Chromosome-based x-axis with color bar (supports naming like 1A, 1B, 1D, etc.)
  - Significance threshold line
  - Optional SNP density coloring

- `plot_pca()`: PCA visualization with 2D/3D support
  - Set `plot_3d = TRUE` for 3D scatter plot (requires scatterplot3d)
  - Group coloring with custom colors
  - Variance explained labels on axes
  - Optional point labels

- `plot_dendrogram()`: Population structure dendrogram
  - Circular (`type = "circular"`) or rectangular layout
  - Group-colored branches and labels
  - Multiple clustering methods (ward.D2, complete, etc.)
  - Requires ape package for circular layout

- `plot_phenotypes()`: Notched boxplots for trait comparison
  - Multi-trait comparison with custom colors
  - Red mean points overlay
  - Notched confidence intervals
  - Publication-ready styling

- `plot_relationship_matrix()`: Improved kinship heatmap
  - Hierarchical clustering option
  - Blue-white-red color scale

### Usage Examples

```r
# Multi-trait Manhattan plot
gwas_list <- list(Yield = gwas1, Height = gwas2, Flowering = gwas3)
plot_manhattan(gwas_list, chr_info = map_data)

# 3D PCA
plot_pca(gdata, groups = populations, plot_3d = TRUE)

# Circular dendrogram
K <- calc_relationship_matrices(gdata, matrices = "A")
plot_dendrogram(gdata, K = K$A, type = "circular", groups = subpopulations)

# Phenotype boxplots
plot_phenotypes(pheno_data, traits = c("Yield", "Height", "Flowering"),
                notch = TRUE, show_mean = TRUE)
```

## New Feature: QTL-Aware GBLUP

Implemented major/minor genetic effect separation. This allows fitting known QTLs as fixed effects while using a reduced kinship matrix for the polygenic background.

### New Functions

- `build_qtl_aware_matrices()`: Build reduced K matrix excluding QTL markers, extract QTL genotypes for fixed effects
- `fit_gblup_qtl()`: Fit GBLUP with QTLs as fixed effects
- `cv_gblup_qtl()`: Cross-validation for QTL-aware models
- `fit_gblup_auto_qtl()`: Automatic pipeline: GWAS → QTL detection → model fitting
- `get_qtl_effects()`: Extract estimated QTL effects from model
- `get_gebv_components()`: Get breakdown of GEBV into polygenic vs QTL components

## New Feature: Strict Cross-Validation (No Data Leakage)

Cross-validation now recomputes the K matrix for each fold using only training genotypes by default. This prevents the subtle data leakage that occurs when test individuals' genotypes contribute to the K matrix.

### New Parameter

- `recompute_K`: Added to `cv_gblup()` function
  - `TRUE` (default): K matrix recomputed per fold from training genotypes only. No data leakage.
  - `FALSE`: K matrix computed once from all individuals. Faster but allows minor genotype leakage.

### Usage Example

```r
# Strict CV (default) - no data leakage
cv_strict <- cv_gblup(gdata, effects = "A", scheme = "CV1", 
                       n_folds = 5, recompute_K = TRUE)

# Fast CV - uses pre-computed K (minor leakage)
K <- calc_relationship_matrices(gdata, matrices = "A")
cv_fast <- cv_gblup(gdata, K, effects = "A", scheme = "CV1", 
                     n_folds = 5, recompute_K = FALSE)
```

### Technical Details

When `recompute_K = TRUE`:
1. Training genotypes (M_train) are extracted for each fold
2. K_train is computed from M_train only using VanRaden formula
3. Model is fitted using K_train
4. Test predictions use proper BLUP equation: û_test = K[test,train] × K[train,train]⁻¹ × û_train
5. K[test,train] is computed from test and training genotypes

This is the statistically correct approach that prevents any information from test individuals leaking into the training process.

## Other Improvements

- Added memory check for Kronecker products in `build_GE_matrix()` 
- Documented MME constants with rationale
- Fixed version display to use dynamic `packageVersion()`
- Added `@importFrom stats cov` for proper namespace imports

### Usage Example (QTL-Aware)

```r
# Method 1: Manual QTL specification
qtl_mats <- build_qtl_aware_matrices(M, qtl_markers = c("SNP_100", "SNP_500"))
model <- fit_gblup_qtl(gdata, Ka = qtl_mats$Ka, Xa = qtl_mats$Xa)

# Method 2: Automatic from GWAS
result <- fit_gblup_auto_qtl(gdata, gwas_threshold = 0.05, 
                              gwas_method = "bonferroni")
print(result$qtl_effects)

# Access QTL effects and GEBV components
get_qtl_effects(model)
get_gebv_components(model)
```

### Key Concepts

- **Ka (reduced K)**: Relationship matrix computed from non-QTL markers
- **Xa (QTL genotypes)**: QTL markers fitted as fixed effects
- **KaeE1**: Per-QTL GxE matrices for environment-specific effects
- **KaeE2**: Residual GxE from non-QTL markers

### When to Use

- Known major genes (disease resistance, quality traits)
- When GWAS detects significant associations
- Oligogenic traits with some large-effect loci
- Want to separate "explained" vs "unexplained" genetic variance

---

# CompGEBLUP 2.0.0

## Major Improvements

### C++ Integration (Previously Unused)
- **Integrated all C++ functions**: The previous version had C++ code that was exported but never called. Now all matrix operations use fast C++ implementations by default.
- **New C++ functions**:
  - `calc_A_matrix_cpp()`: Fast VanRaden A matrix calculation using RcppArmadillo
  - `calc_D_matrix_cpp()`: Fast Vitezica dominance matrix calculation
  - `calc_epistatic_matrix_cpp()`: Fast Hadamard product for epistatic matrices
  - `emma_reml_cpp()`: EMMA-style efficient REML with O(n) per-iteration complexity
  - `gwas_emmax_cpp()`: EMMAX-style fast GWAS marker testing
  - `eigen_decomp_cpp()`: Fast eigendecomposition
  - `create_Z_sparse_cpp()`: Sparse incidence matrix creation
  - `compute_ZKZt_cpp()`: Efficient ZKZ' computation
  - `solve_mme_cpp()`: Direct MME solver
  - `make_pd_cpp()`: Fast positive definite matrix adjustment
  - `is_positive_definite_cpp()`: Fast PD check

### Statistical Methods

#### EMMA-style Efficient Algorithm
- Implemented spectral decomposition approach for single random effect models
- Reduces per-iteration complexity from O(n³) to O(n)
- Uses golden section search for λ optimization
- Automatically selected when appropriate (single random effect)

#### Fixed Formula Inconsistencies
- Unified VanRaden (2008) A matrix formula between R and C++
- Unified Vitezica (2013) D matrix formula
- Consistent denominator calculations across implementations

#### Proper BLUP Prediction
- Fixed GEBV prediction for untrained individuals using proper BLUP equations:
  ```
  u_new = K_new,train * K_train,train^{-1} * u_train
  ```
- New `predict_blup_new()` function for proper out-of-sample prediction
- Removed heuristic K-matrix interpolation

### Performance Optimizations

#### Sparse Matrix Support
- Added optional sparse matrix support for design matrices (Z)
- New parameter `use.sparse` in `fit_gblup()` and `build_design_matrices()`
- Uses `Matrix::sparseMatrix` for memory-efficient storage
- Particularly beneficial for datasets with n > 1000

#### Memory Efficiency
- Reduced memory footprint for large datasets
- Improved eigendecomposition caching
- Vectorized operations where possible

#### OpenMP Parallelization
- GWAS marker testing supports OpenMP parallel processing
- Enabled in Makevars with `PKG_CXXFLAGS = $(SHLIB_OPENMP_CXXFLAGS)`

### Code Quality Improvements

#### Replaced Magic Numbers with Constants
- Created `.MME_CONSTANTS` list with documented values:
  - `MIN_VARIANCE = 1e-10`
  - `MIN_EIGENVALUE = 1e-10`
  - `DEFAULT_RIDGE = 1e-6`
  - `MIN_OBSERVATIONS = 5`
  - `STEP_HALVING_INIT = 0.9`
  - `STEP_HALVING_MIN = 0.01`

#### Improved Error Handling
- Consistent use of `stop()` for errors, `warning()` for warnings
- More informative error messages
- Better handling of edge cases (singular matrices, zero variance)

#### Fixed Division by Zero
- Proper handling of near-zero denominators in A and D matrix calculations
- Warning issued before falling back to safe value

### API Improvements

#### New Parameters
- `use.cpp`: Control whether to use C++ implementations (default TRUE)
- `use.emma`: Control whether to use EMMA algorithm (default TRUE)
- `use.sparse`: Control whether to use sparse matrices (default FALSE)

#### Deprecated Parameters
- `ridge` parameter in `make_positive_definite()` deprecated in favor of `tol`

#### New Functions
- `predict_blup_new()`: Proper BLUP prediction for new individuals
- `compute_pev()`: Prediction error variance calculation
- `compute_reliability()`: Reliability of GEBVs
- `is_positive_definite()`: Check if matrix is PD

### Documentation

- Updated all function documentation
- Added examples for C++ functions
- Improved vignette with performance comparisons

## Bug Fixes

- Fixed filter function name collision (`filter_maf` → `filter_markers_maf`)
- Fixed character/factor conversion issues throughout
- Fixed GEBV extraction for environment-specific predictions
- Fixed lambda calculation threshold (increased to 100 markers)

## Breaking Changes

- C++ function names changed:
  - `calc_A_cpp` → `calc_A_matrix_cpp`
  - `calc_D_cpp` → `calc_D_matrix_cpp`
  - `hadamard_product` → `calc_epistatic_matrix_cpp`
- `filter_maf()` renamed to `filter_markers_maf()` (internal function)

## Performance Benchmarks

Approximate speedup factors (vs v1.0):

| Operation | n=500 | n=1000 | n=2000 |
|-----------|-------|--------|--------|
| A matrix  | 2x    | 3x     | 5x     |
| D matrix  | 2x    | 3x     | 5x     |
| REML (single RE) | 5x | 10x | 20x |
| GWAS (1000 markers) | 3x | 5x | 8x |

---

# CompGEBLUP 1.0.0

Initial release with:
- GBLUP model fitting
- Cross-validation (CV0, CV1, CV2)
- GWAS
- Multi-trait analysis
- 5 pre-defined models
