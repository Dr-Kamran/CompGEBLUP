# CompGEBLUP

**Comprehensive Genomic Prediction with Multi-Trait GBLUP**

[![R](https://img.shields.io/badge/R-%E2%89%A54.0.0-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-GPL%20(%E2%89%A53)-brightgreen.svg)](https://www.gnu.org/licenses/gpl-3.0)

A complete toolkit for genomic selection including additive, dominance, epistatic effects, G×E interactions, multi-trait selection, cross-validation, and GWAS. Features an independent REML-based mixed model engine with EMMA-style efficient algorithms and fast C++ matrix operations using RcppArmadillo.

## ✨ Key Features

- **Fast PLINK BED Loading**: C++ backend for ~1000× speedup on large datasets
- **Relationship Matrices**: Additive (VanRaden 2008), Dominance (Vitezica 2013), Epistatic (A×A, A×D, D×D)
- **Variance Component Estimation**: AI-REML for multiple effects, EMMA for single effect
- **Cross-Validation**: k-fold CV with proper data leakage prevention
- **G×E Interactions**: Multi-environment modeling
- **GWAS**: EMMAX-style mixed model association
- **QTL-aware GBLUP**: Separate major QTL from polygenic background
- **Multi-trait Selection**: Correlated trait prediction

## 🚀 Quick Start

### Installation

```r
# From GitHub
devtools::install_github("Dr-Kamran/CompGEBLUP")

# From local repository
devtools::install()  # Run from package root
```

### Basic Usage

```r
library(CompGEBLUP)

# Load genotypes (NEW: Fast PLINK loader!)
plink_data <- load_plink_bed("my_data.bed", coding = "additive")
M <- plink_data$geno  # Matrix ready for analysis

# Create dataset
pheno <- data.frame(GID = rownames(M), ENV = "E1", trait = your_phenotypes)
gdata <- new("GBLUPData", genotypes = M, phenotypes = pheno)

# Build relationship matrix
A <- build_A_matrix(M, method = "VanRaden")

# Fit GBLUP model
model <- fit_gblup(gdata, K_matrices = list(A = A), effects = "A")

# Get heritability
h2 <- heritability(model)
print(h2)

# Cross-validate
cv <- cv_gblup(gdata, list(A = A), "A", n_folds = 5)
print(cv@metrics)
```

## 📚 Documentation

- **[USAGE_GUIDE.md](USAGE_GUIDE.md)**: Complete installation and usage guide
- **[QUICK_REFERENCE.md](QUICK_REFERENCE.md)**: Quick reference for common tasks
- **[G2F Analysis Guide](inst/examples/README.md)**: Detailed guide for hybrid maize analysis
- **[G2F Example Script](inst/examples/G2F_analysis.R)**: Complete working example

### Function Help

```r
?load_plink_bed     # Fast PLINK BED file loading
?build_A_matrix     # Additive relationship matrix
?build_D_matrix     # Dominance relationship matrix
?fit_gblup          # Fit GBLUP model
?cv_gblup           # Cross-validation
?gwas               # Genome-wide association
```

## 🆕 What's New: Fast PLINK Loader

The package now includes a C++ backend for loading PLINK binary files:

```r
# Load large datasets quickly (1000 ind × 500k SNPs in seconds)
data <- load_plink_bed("genotypes.bed", 
                       coding = "additive",    # For CompGEBLUP
                       impute_missing = TRUE,
                       verbose = TRUE)

M <- data$geno  # Genotype matrix
fam <- data$fam # Individual metadata
bim <- data$bim # SNP metadata
```

**Performance**: ~1000× faster than pure R implementations for large datasets.

## 📊 Example Workflows

### 1. Basic Genomic Prediction

```r
# Load and prepare data
M <- load_plink_bed("data.bed", coding = "additive")$geno
gdata <- new("GBLUPData", genotypes = M, phenotypes = pheno_df)

# Build matrix and fit model
A <- build_A_matrix(M)
model <- fit_gblup(gdata, list(A = A), "A")

# Results
summary(model)
heritability(model)
```

### 2. Model Comparison (A vs A+D)

```r
# Build matrices
A <- build_A_matrix(M)
D <- build_D_matrix(M)

# Fit models
model_A <- fit_gblup(gdata, list(A = A), "A")
model_AD <- fit_gblup(gdata, list(A = A, D = D), c("A", "D"))

# Cross-validate
cv_A <- cv_gblup(gdata, list(A = A), "A", n_folds = 5, n_reps = 3)
cv_AD <- cv_gblup(gdata, list(A = A, D = D), c("A", "D"), n_folds = 5, n_reps = 3)

# Compare
data.frame(
  Model = c("A", "A+D"),
  Accuracy = c(mean(cv_A@metrics$predictive_ability),
               mean(cv_AD@metrics$predictive_ability))
)
```

### 3. Hybrid Maize Analysis (G2F)

See **[inst/examples/G2F_analysis.R](inst/examples/G2F_analysis.R)** for a complete workflow including:
- Inbred parent genotype loading
- Hybrid genotype construction
- Multi-environment phenotype handling
- Outlier filtering
- Variance component estimation
- Cross-validation
- Model comparison (A, A+D, A+D+AA)

## 🔬 Methods Implemented

### Relationship Matrices
- **Additive (A)**: VanRaden (2008), UAR, UARadj
- **Dominance (D)**: Vitezica et al. (2013)
- **Epistatic (AA, AD, DD)**: Hadamard products

### Variance Component Estimation
- **EMMA**: Kang et al. (2008) - Efficient for single random effect
- **AI-REML**: Gilmour et al. (1995) - Robust for multiple effects

### Cross-Validation Schemes
- **CV1**: Random k-fold
- **CV2**: Leave-one-environment-out
- **CV0**: Within-environment k-fold

## 📦 Input Data Formats

### PLINK Files
- `.bed`: Binary genotypes (SNP-major mode)
- `.bim`: SNP information (6 columns: CHR, SNP, CM, BP, A1, A2)
- `.fam`: Sample information (6 columns: FID, IID, PID, MID, Sex, Pheno)

### Phenotype Data
```r
pheno_df <- data.frame(
  GID = c("Ind1", "Ind1", "Ind2", "Ind2"),
  ENV = c("E1", "E2", "E1", "E2"),
  trait = c(100.5, 98.3, 105.2, 102.8)
)
```

### Hybrid Phenotype Data (for G2F analysis)
```csv
hybrid_id,environment,trait,parent1,parent2
H001,Location1,12.5,P001,P002
H001,Location2,11.8,P001,P002
```

## 🛠️ Requirements

- R ≥ 4.0.0
- Rcpp ≥ 1.0.0
- RcppArmadillo
- Matrix
- Rtools (Windows) or build tools (Mac/Linux) for C++ compilation

## 📈 Performance

| Operation | Dataset Size | Time |
|-----------|--------------|------|
| PLINK BED loading | 1000 ind × 500k SNPs | ~5 seconds |
| A matrix (C++) | 1000 × 10k SNPs | ~2 seconds |
| EMMA REML (1 effect) | 1000 individuals | ~5 seconds |
| AI-REML (2 effects) | 1000 individuals | ~30 seconds |

## 🐛 Troubleshooting

| Issue | Solution |
|-------|----------|
| Package won't compile | Install Rcpp, RcppArmadillo; check Rtools |
| BED file error | Verify with PLINK: `plink --bfile data --make-bed` |
| Memory issues | Filter SNPs by MAF, reduce cross-validation folds |
| Convergence warnings | Increase `nIters` to 500, adjust `ridge` parameter |

See **[USAGE_GUIDE.md](USAGE_GUIDE.md)** for detailed troubleshooting.

## 📖 References

1. **VanRaden PM (2008)** Efficient methods to compute genomic predictions. *Journal of Dairy Science* 91:4414-4423.

2. **Vitezica ZG et al. (2013)** On the additive and dominant variance and covariance of individuals within the genomic selection scope. *Genetics* 195:1223-1230.

3. **Kang HM et al. (2008)** Efficient control of population structure in model organism association mapping. *Genetics* 178:1709-1723.

4. **Gilmour AR et al. (1995)** Average information REML: An efficient algorithm for variance parameter estimation in linear mixed models. *Biometrics* 51:1440-1450.

## 👥 Author

**Muhammad Kamran**  
Email: kamran@zju.edu.cn  
ORCID: [0009-0006-1763-7536](https://orcid.org/0009-0006-1763-7536)

## 📄 License

GPL (≥ 3)

## 🤝 Contributing

Contributions are welcome! Please:
1. Fork the repository
2. Create a feature branch
3. Make your changes with tests
4. Submit a pull request

## 🔗 Links

- **GitHub Repository**: https://github.com/Dr-Kamran/CompGEBLUP
- **Issues**: https://github.com/Dr-Kamran/CompGEBLUP/issues
- **Documentation**: See [USAGE_GUIDE.md](USAGE_GUIDE.md)

## 📝 Citation

If you use CompGEBLUP in your research, please cite:

```bibtex
@software{CompGEBLUP,
  author = {Kamran, Muhammad},
  title = {CompGEBLUP: Comprehensive Genomic Prediction with Multi-Trait GBLUP},
  year = {2024},
  url = {https://github.com/Dr-Kamran/CompGEBLUP}
}
```

---

**Getting Started**: See [USAGE_GUIDE.md](USAGE_GUIDE.md) for installation and usage instructions.

**Quick Reference**: See [QUICK_REFERENCE.md](QUICK_REFERENCE.md) for command cheat sheet.

**G2F Example**: See [inst/examples/README.md](inst/examples/README.md) for hybrid maize analysis guide.
