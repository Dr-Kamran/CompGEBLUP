# G2F CompGEBLUP Analysis Script

This directory contains a corrected and optimized analysis script for Genomes to Fields (G2F) hybrid maize data using the CompGEBLUP package.

## Overview

The `G2F_analysis.R` script demonstrates best practices for genomic prediction in hybrid crops, addressing multiple statistical and computational issues found in standalone scripts.

## Key Features

### 1. **Fast PLINK BED File Loading**
- Uses C++ backend (`decode_bed_cpp`) for efficient genotype decoding
- Handles large datasets (1000+ individuals × 500,000+ SNPs) in seconds
- Properly validates BED file format and magic numbers

### 2. **Correct Statistical Methods**
- **Additive matrix (A)**: VanRaden (2008) method with C++ acceleration
- **Dominance matrix (D)**: Vitezica et al. (2013) method (NOT Su et al. 2012)
- **Epistatic matrix (AA)**: Hadamard product A#A
- **Hybrid D matrix**: Correctly uses D matrix entries for diagonal

### 3. **Robust Variance Component Estimation**
- AI-REML algorithm for multiple random effects (A+D, A+D+AA models)
- EMMA algorithm for single random effect (A model) - very fast
- Proper convergence checking with log-likelihood monitoring
- Up to 500 iterations with step-halving when needed

### 4. **Cross-Validation**
- k-fold CV with recomputed K matrices per fold (prevents data leakage)
- Multiple replications for stable estimates
- Reports predictive ability (correlation) and MSE

### 5. **Data Quality Controls**
- Phenotype outlier filtering (±4 SD within environment)
- Multi-environment data handling (no information loss)
- Missing data imputation

### 6. **Reproducibility**
- Random seeds set for fold assignment
- Memory management with explicit `gc()` calls
- Clear comments explaining each step

## Usage

### Prerequisites

Install the CompGEBLUP package:
```r
# From GitHub
# devtools::install_github("Dr-Kamran/CompGEBLUP")

library(CompGEBLUP)
```

### Input Data Requirements

#### 1. Genotype Data (PLINK format)
Three files with the same prefix:
- `*.bed` - Binary genotype file
- `*.bim` - Variant information (6 columns: chr, rsID, genetic distance, position, allele1, allele2)
- `*.fam` - Sample information (6 columns: FID, IID, paternal ID, maternal ID, sex, phenotype)

The BED file must be in SNP-major mode (standard for PLINK 1.9).

#### 2. Phenotype Data (CSV format)
Required columns:
- `hybrid_id` - Unique hybrid identifier
- `environment` - Environment/location code
- `trait` - Phenotype value (e.g., grain yield)
- `parent1` - Inbred parent 1 ID (must match `*.fam` IID)
- `parent2` - Inbred parent 2 ID (must match `*.fam` IID)

Example:
```csv
hybrid_id,environment,trait,parent1,parent2
H001,IA_2019,12.5,P001,P002
H001,IL_2019,11.8,P001,P002
H002,IA_2019,13.2,P003,P004
```

### Running the Analysis

1. **Edit the script to set your data paths:**
```r
# Line ~50
plink_prefix <- "/path/to/your/data/g2f_genotypes"

# Line ~80
pheno_file <- "/path/to/your/data/phenotypes.csv"
```

2. **Run the script:**
```r
source("inst/examples/G2F_analysis.R")
```

The script will:
1. Load and validate genotype data
2. Load and filter phenotype data
3. Construct hybrid genotype matrix from parent genotypes
4. Build relationship matrices (A, D, AA)
5. Estimate variance components for 3 models
6. Perform cross-validation
7. Report model comparison

### Expected Runtime

For a typical G2F dataset:
- 500 inbred lines × 300,000 SNPs: ~2-5 minutes loading
- 2,000 hybrids × 5 environments: ~10-30 minutes for full analysis
- Cross-validation (5-fold × 3 reps): ~20-60 minutes per model

Total: ~1-2 hours for complete analysis with cross-validation.

## Output

The script produces:

### Variance Components
```
Model 1: Additive (A)
  σ²_A = 12.5
  σ²_e = 8.3
  h² = 0.60

Model 2: Additive + Dominance (A+D)
  σ²_A = 10.2
  σ²_D = 2.8
  σ²_e = 8.1
  h²_narrow = 0.48
  h²_broad = 0.62
```

### Cross-Validation Results
```
Model A Cross-Validation Results:
  Mean Predictive Ability: 0.523 ± 0.042
  Mean MSE: 6.82 ± 0.54

Model A+D Cross-Validation Results:
  Mean Predictive Ability: 0.547 ± 0.038
  Mean MSE: 6.45 ± 0.51
```

### Model Comparison Table
```
  Model  Heritability  Predictive_Ability   MSE
      A         0.600              0.523  6.82
    A+D         0.620              0.547  6.45
 A+D+AA         0.635              0.551  6.38
```

## Key Improvements Over Standalone Scripts

| Issue | Old Approach | New Approach | Speedup |
|-------|--------------|--------------|---------|
| **BED parsing** | Pure R loops | C++ `decode_bed_cpp` | 1000× |
| **D matrix** | Su et al. (2012) | Vitezica et al. (2013) | Correct method |
| **Hybrid D diagonal** | `1 - A[p1,p2]` (wrong) | Computed correctly by `build_D_matrix()` | Correct formula |
| **Variance estimation** | Naive EM (100 iter) | AI-REML (500 iter) | More robust |
| **Convergence** | No checks | Log-likelihood + step-halving | Reliable |
| **Cross-validation** | None | k-fold with K recomputation | Unbiased |
| **Outliers** | None | ±4 SD filtering | Better QC |
| **Multi-environment** | Averaged | Individual records | No G×E loss |

## Customization

### Adjusting Cross-Validation
```r
# Change number of folds
cv_A <- cv_gblup(..., n_folds = 10, ...)

# Change number of replications
cv_A <- cv_gblup(..., n_reps = 5, ...)

# Use different CV scheme
cv_A <- cv_gblup(..., scheme = "CV2", ...)  # Leave-one-environment-out
```

### Adding More Models
```r
# Fit A+D+AD model
K_ADAD <- list(A = A, D = D, AD = build_E_matrix(M_hybrids, type = "A#D"))
model_ADAD <- fit_gblup(gdata, K_ADAD, effects = c("A", "D", "AD"))
```

### Prediction for New Hybrids
See Section 6 in the script for prediction workflow.

## References

1. **VanRaden PM (2008)** Efficient methods to compute genomic predictions. *Journal of Dairy Science* 91:4414-4423.

2. **Vitezica ZG et al. (2013)** On the additive and dominant variance and covariance of individuals within the genomic selection scope. *Genetics* 195:1223-1230.

3. **Su G et al. (2012)** Estimating additive and non-additive genetic variances and predicting genetic merits using genome-wide dense single nucleotide polymorphism markers. *PLoS ONE* 7:e45293.

4. **Kang HM et al. (2008)** Efficient control of population structure in model organism association mapping. *Genetics* 178:1709-1723.

5. **Gilmour AR et al. (1995)** Average information REML: An efficient algorithm for variance parameter estimation in linear mixed models. *Biometrics* 51:1440-1450.

## Troubleshooting

### Memory Issues
```r
# Reduce number of SNPs via filtering
M <- M[, maf > 0.05]  # Keep only common variants

# Use faster CV (allows minor data leakage)
cv_A <- cv_gblup(..., recompute_K = FALSE, ...)
```

### Convergence Problems
```r
# Increase ridge parameter
model <- fit_gblup(..., ridge = 0.01, ...)

# Increase iterations
model <- fit_gblup(..., nIters = 1000, ...)
```

### Missing Parent Genotypes
The script automatically removes hybrids whose parent genotypes are not available. Ensure parent IDs in the phenotype file exactly match the IIDs in the `.fam` file.

## Contact

For questions about the CompGEBLUP package:
- GitHub Issues: https://github.com/Dr-Kamran/CompGEBLUP/issues
- Author: Muhammad Kamran (kamran@zju.edu.cn)
