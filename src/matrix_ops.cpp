// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace Rcpp;
using namespace arma;

// Constants for numerical stability (match R implementation)
const double MIN_VARIANCE = 1e-10;
const double DEFAULT_RIDGE = 1e-6;
const double RIDGE_MULTIPLIER = 10.0;

//' Fast Additive Relationship Matrix (VanRaden 2008)
//' 
//' Calculates the genomic relationship matrix using the VanRaden (2008) method.
//' Uses RcppArmadillo for optimized linear algebra operations.
//' 
//' @param G Genotype matrix (individuals x markers) coded as -1, 0, 1
//' @return Additive relationship matrix (n x n)
//' @export
// [[Rcpp::export]]
arma::mat calc_A_matrix_cpp(const arma::mat& G) {
  int p = G.n_cols;
  
  // Calculate allele frequencies (convert -1,0,1 to 0,1,2 then get freq)
  arma::rowvec maf = arma::mean(G + 1.0, 0) / 2.0;
  
  // Center genotypes: W = G - 2*(p - 0.5) = G - (2p - 1)
  arma::mat W = G;
  for (int j = 0; j < p; j++) {
    double center = 2.0 * maf(j) - 1.0;
    W.col(j) -= center;
  }
  
  // Calculate denominator: 2 * sum(p * (1-p))
  double denom = 2.0 * arma::accu(maf % (1.0 - maf));
  
  // Prevent division by zero
  if (denom < 1e-10) {
    Rcpp::warning("Denominator near zero in A matrix calculation. Using 1.0");
    denom = 1.0;
  }
  
  // A = WW' / denom (using efficient matrix multiplication)
  arma::mat A = (W * W.t()) / denom;
  
  // Ensure symmetry
  A = (A + A.t()) / 2.0;
  
  return A;
}

//' Fast Dominance Relationship Matrix (Vitezica 2013)
//' 
//' Calculates the dominance relationship matrix using Vitezica et al. (2013) method.
//' 
//' @param G Genotype matrix (individuals x markers) coded as -1, 0, 1
//' @return Dominance relationship matrix (n x n)
//' @export
// [[Rcpp::export]]
arma::mat calc_D_matrix_cpp(const arma::mat& G) {
  int n = G.n_rows;
  int p = G.n_cols;
  
  // Calculate allele frequencies
  arma::rowvec maf = arma::mean(G + 1.0, 0) / 2.0;
  
  // Create heterozygote indicator matrix (1 if G == 0, else 0)
  arma::mat H = arma::zeros<arma::mat>(n, p);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < p; j++) {
      if (std::abs(G(i, j)) < 0.5) {  // G == 0 (heterozygote)
        H(i, j) = 1.0;
      }
    }
  }
  
  // Expected heterozygosity: 2*p*(1-p)
  arma::rowvec het_exp = 2.0 * maf % (1.0 - maf);
  
  // Center: W = H - het_exp
  arma::mat W = H;
  for (int j = 0; j < p; j++) {
    W.col(j) -= het_exp(j);
  }
  
  // Denominator: sum((2*p*(1-p))^2) = sum(het_exp^2)
  double denom = arma::accu(arma::square(het_exp));
  
  if (denom < 1e-10) {
    Rcpp::warning("Denominator near zero in D matrix calculation. Using 1.0");
    denom = 1.0;
  }
  
  // D = WW' / denom
  arma::mat D = (W * W.t()) / denom;
  
  // Ensure symmetry
  D = (D + D.t()) / 2.0;
  
  return D;
}

//' Fast Epistatic Relationship Matrix (Hadamard Product)
//' 
//' Calculates epistatic relationship matrix as element-wise product of two matrices.
//' 
//' @param A First relationship matrix
//' @param B Second relationship matrix
//' @return Epistatic relationship matrix (Hadamard product)
//' @export
// [[Rcpp::export]]
arma::mat calc_epistatic_matrix_cpp(const arma::mat& A, const arma::mat& B) {
  if (A.n_rows != B.n_rows || A.n_cols != B.n_cols) {
    Rcpp::stop("Matrices must have same dimensions");
  }
  
  // Element-wise multiplication (Hadamard product)
  arma::mat E = A % B;
  
  // Ensure symmetry
  E = (E + E.t()) / 2.0;
  
  return E;
}

//' Eigendecomposition for EMMA-style efficient computation
//' 
//' Performs spectral decomposition of the relationship matrix for 
//' efficient variance component estimation.
//' 
//' @param K Relationship matrix (symmetric positive semi-definite)
//' @return List with eigenvectors (U) and eigenvalues (d)
//' @export
// [[Rcpp::export]]
Rcpp::List eigen_decomp_cpp(const arma::mat& K) {
  // Ensure input is symmetric
  arma::mat K_sym = (K + K.t()) / 2.0;
  
  arma::vec eigval;
  arma::mat eigvec;
  
  // Symmetric eigendecomposition
  arma::eig_sym(eigval, eigvec, K_sym);
  
  // Sort in descending order (arma gives ascending)
  arma::uvec idx = arma::sort_index(eigval, "descend");
  eigval = eigval(idx);
  eigvec = eigvec.cols(idx);
  
  // Floor small eigenvalues to prevent numerical issues
  int n_floored = 0;
  for (size_t i = 0; i < eigval.n_elem; i++) {
    if (eigval(i) < 1e-10) {
      eigval(i) = 1e-10;
      n_floored++;
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("vectors") = eigvec,
    Rcpp::Named("values") = eigval,
    Rcpp::Named("n_floored") = n_floored
  );
}

//' Fast REML log-likelihood evaluation (EMMA-style)
//' 
//' Computes the REML log-likelihood given transformed data and lambda.
//' Uses pre-computed eigendecomposition for O(n) complexity per evaluation.
//' 
//' @param eta Transformed response (U'y where U is eigenvector matrix)
//' @param xi Transformed fixed effects (U'X)
//' @param d Eigenvalues of K matrix
//' @param lambda Ratio of residual to genetic variance (sigma2_e / sigma2_g)
//' @return REML log-likelihood value
//' @export
// [[Rcpp::export]]
double emma_reml_loglik_cpp(const arma::vec& eta, const arma::mat& xi,
                            const arma::vec& d, double lambda) {
  // Diagonal of V in transformed space: d + lambda
  arma::vec D_lambda = d + lambda;
  
  // Inverse weights
  arma::vec D_inv = 1.0 / D_lambda;
  
  // Weighted X'X and X'y
  arma::mat xiT_D_xi = xi.t() * (xi.each_col() % D_inv);
  arma::vec xiT_D_eta = xi.t() * (eta % D_inv);
  
  // Solve for beta
  arma::vec beta;
  bool solved = arma::solve(beta, xiT_D_xi, xiT_D_eta);
  if (!solved) {
    // Add small ridge if singular
    xiT_D_xi.diag() += 1e-8;
    arma::solve(beta, xiT_D_xi, xiT_D_eta);
  }
  
  // Residuals
  arma::vec resid = eta - xi * beta;
  
  // Weighted sum of squares
  double rss = arma::accu(arma::square(resid) % D_inv);
  
  // Log determinants
  double log_det_V = arma::accu(arma::log(D_lambda));
  double log_det_XVX;
  double sign;
  arma::log_det(log_det_XVX, sign, xiT_D_xi);
  
  // REML log-likelihood
  double loglik = -0.5 * (log_det_V + log_det_XVX + rss);
  
  // Return negative infinity for invalid results
  if (!std::isfinite(loglik)) {
    return -std::numeric_limits<double>::infinity();
  }
  
  return loglik;
}

//' EMMA-style variance component estimation using Brent's method
//' 
//' Efficiently estimates variance components using spectral decomposition
//' and 1D optimization over lambda = sigma2_e / sigma2_g.
//' 
//' @param y Response vector
//' @param X Fixed effects design matrix
//' @param K Relationship matrix
//' @param tol Convergence tolerance
//' @param max_iter Maximum iterations for optimization
//' @return List with sigma2_g, sigma2_e, and other parameters
//' @export
// [[Rcpp::export]]
Rcpp::List emma_reml_cpp(const arma::vec& y, const arma::mat& X,
                         const arma::mat& K, double tol = 1e-6, 
                         int max_iter = 100) {
  int n = y.n_elem;
  int p = X.n_cols;
  
  // Ensure K is symmetric (handle floating point errors)
  arma::mat K_sym = (K + K.t()) / 2.0;
  
  // Eigendecomposition of K
  arma::vec d;
  arma::mat U;
  arma::eig_sym(d, U, K_sym);
  
  // Floor small eigenvalues
  for (size_t i = 0; i < d.n_elem; i++) {
    if (d(i) < 1e-10) d(i) = 1e-10;
  }
  
  // Transform data: eta = U'y, xi = U'X
  arma::vec eta = U.t() * y;
  arma::mat xi = U.t() * X;
  
  // Grid search + Brent's method for lambda optimization
  // Start with grid search
  arma::vec lambda_grid = arma::logspace(-5, 5, 20);
  double best_lambda = 1.0;
  double best_ll = -arma::datum::inf;
  
  for (size_t i = 0; i < lambda_grid.n_elem; i++) {
    double ll = emma_reml_loglik_cpp(eta, xi, d, lambda_grid(i));
    if (ll > best_ll) {
      best_ll = ll;
      best_lambda = lambda_grid(i);
    }
  }
  
  // Refine with golden section search
  double a = best_lambda * 0.1;
  double b = best_lambda * 10.0;
  double gr = (std::sqrt(5.0) - 1.0) / 2.0;  // golden ratio
  
  double c = b - gr * (b - a);
  double dd = a + gr * (b - a);
  
  for (int iter = 0; iter < max_iter; iter++) {
    double fc = emma_reml_loglik_cpp(eta, xi, d, c);
    double fd = emma_reml_loglik_cpp(eta, xi, d, dd);
    
    if (fc > fd) {
      b = dd;
      dd = c;
      c = b - gr * (b - a);
    } else {
      a = c;
      c = dd;
      dd = a + gr * (b - a);
    }
    
    if (std::abs(b - a) < tol * (std::abs(c) + std::abs(dd))) {
      break;
    }
  }
  
  double lambda_opt = (a + b) / 2.0;
  
  // Estimate variance components
  arma::vec D_inv = 1.0 / (d + lambda_opt);
  
  // Solve for beta
  arma::mat xiT_D_xi = xi.t() * (xi.each_col() % D_inv);
  arma::vec xiT_D_eta = xi.t() * (eta % D_inv);
  arma::vec beta;
  arma::solve(beta, xiT_D_xi, xiT_D_eta);
  
  // Residuals and variance estimates
  arma::vec resid = eta - xi * beta;
  double rss = arma::accu(arma::square(resid) % D_inv);
  
  // REML estimates
  double sigma2_e = rss / (n - p);
  double sigma2_g = sigma2_e / lambda_opt;
  
  // Heritability
  double h2 = sigma2_g / (sigma2_g + sigma2_e);
  
  // Final log-likelihood
  double loglik = emma_reml_loglik_cpp(eta, xi, d, lambda_opt);
  
  return Rcpp::List::create(
    Rcpp::Named("sigma2_g") = sigma2_g,
    Rcpp::Named("sigma2_e") = sigma2_e,
    Rcpp::Named("lambda") = lambda_opt,
    Rcpp::Named("h2") = h2,
    Rcpp::Named("beta") = beta,
    Rcpp::Named("loglik") = loglik,
    Rcpp::Named("U") = U,
    Rcpp::Named("d") = d
  );
}

//' Fast GWAS marker testing with pre-computed decomposition
//' 
//' Tests all markers efficiently using pre-computed eigendecomposition.
//' This is the core of EMMAX/GEMMA-style fast GWAS.
//' 
//' @param y Response vector
//' @param X Base fixed effects design matrix (usually just intercept)
//' @param M Marker matrix (n x m)
//' @param U Eigenvector matrix from K decomposition
//' @param d Eigenvalues from K decomposition
//' @param lambda Ratio sigma2_e/sigma2_g (from null model)
//' @return Matrix with columns: beta, se, t_stat, p_value
//' @export
// [[Rcpp::export]]
arma::mat gwas_emmax_cpp(const arma::vec& y, const arma::mat& X,
                         const arma::mat& M, const arma::mat& U,
                         const arma::vec& d, double lambda) {
  int n = y.n_elem;
  int m = M.n_cols;
  int p = X.n_cols;
  
  // Transform response and fixed effects
  arma::vec eta = U.t() * y;
  arma::mat xi = U.t() * X;
  
  // Inverse weights
  arma::vec D_inv = 1.0 / (d + lambda);
  
  // Results matrix
  arma::mat results(m, 4);  // beta, se, t_stat, p_value
  
  // Degrees of freedom
  double df = n - p - 1;
  if (df < 1) df = 1;
  
  // Test each marker
  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (int j = 0; j < m; j++) {
    // Transform marker
    arma::vec marker_trans = U.t() * M.col(j);
    
    // Full design matrix with marker
    arma::mat xi_full = arma::join_rows(xi, marker_trans);
    
    // Weighted least squares
    arma::mat xiT_D_xi = xi_full.t() * (xi_full.each_col() % D_inv);
    arma::vec xiT_D_eta = xi_full.t() * (eta % D_inv);
    
    // Solve
    arma::vec beta_full;
    bool solved = arma::solve(beta_full, xiT_D_xi, xiT_D_eta);
    
    if (!solved || !beta_full.is_finite()) {
      results(j, 0) = arma::datum::nan;
      results(j, 1) = arma::datum::nan;
      results(j, 2) = arma::datum::nan;
      results(j, 3) = arma::datum::nan;
      continue;
    }
    
    // Marker effect is last coefficient
    double beta_marker = beta_full(p);
    
    // Compute residual variance for SE
    arma::vec resid = eta - xi_full * beta_full;
    double rss = arma::accu(arma::square(resid) % D_inv);
    double sigma2_hat = rss / df;
    
    // Invert for SE
    arma::mat cov_beta;
    bool inv_success = arma::inv_sympd(cov_beta, xiT_D_xi);
    if (!inv_success) {
      xiT_D_xi.diag() += 1e-8;
      arma::inv_sympd(cov_beta, xiT_D_xi);
    }
    
    double se_marker = std::sqrt(sigma2_hat * cov_beta(p, p));
    
    // T-statistic and p-value
    double t_stat = beta_marker / se_marker;
    double p_val = 2.0 * R::pt(std::abs(t_stat), df, 0, 0);
    
    results(j, 0) = beta_marker;
    results(j, 1) = se_marker;
    results(j, 2) = t_stat;
    results(j, 3) = p_val;
  }
  
  return results;
}

//' Sparse incidence matrix creation
//' 
//' Creates a sparse incidence (Z) matrix mapping observations to levels.
//' 
//' @param factor_idx Integer vector of factor levels (1-indexed)
//' @param n_levels Number of unique levels
//' @return Sparse matrix (dgCMatrix format)
//' @export
// [[Rcpp::export]]
arma::sp_mat create_Z_sparse_cpp(const arma::uvec& factor_idx, int n_levels) {
  int n = factor_idx.n_elem;
  
  // Create sparse matrix
  arma::sp_mat Z(n, n_levels);
  
  for (int i = 0; i < n; i++) {
    int j = factor_idx(i) - 1;  // Convert to 0-indexed
    if (j >= 0 && j < n_levels) {
      Z(i, j) = 1.0;
    }
  }
  
  return Z;
}

//' Efficient computation of ZKZ'
//' 
//' Computes Z * K * Z' efficiently for variance component models.
//' 
//' @param Z Incidence matrix (can be sparse)
//' @param K Relationship matrix (dense)
//' @return ZKZ' matrix
//' @export
// [[Rcpp::export]]
arma::mat compute_ZKZt_cpp(const arma::mat& Z, const arma::mat& K) {
  // Z * K
  arma::mat ZK = Z * K;
  
  // (ZK) * Z'
  arma::mat ZKZt = ZK * Z.t();
  
  return ZKZt;
}

//' Solve mixed model equations efficiently
//' 
//' Solves the Henderson's mixed model equations using efficient algorithms.
//' 
//' @param y Response vector
//' @param X Fixed effects design matrix
//' @param Z Random effects design matrix
//' @param K Relationship matrix
//' @param sigma2_g Genetic variance
//' @param sigma2_e Residual variance
//' @return List with fixed effects (b) and random effects (u)
//' @export
// [[Rcpp::export]]
Rcpp::List solve_mme_cpp(const arma::vec& y, const arma::mat& X,
                         const arma::mat& Z, const arma::mat& K,
                         double sigma2_g, double sigma2_e) {
  int n = y.n_elem;
  
  // Build V = sigma2_g * ZKZ' + sigma2_e * I
  arma::mat ZKZt = Z * K * Z.t();
  arma::mat V = sigma2_g * ZKZt + sigma2_e * arma::eye(n, n);
  
  // Invert V (use Cholesky if PD)
  arma::mat V_inv;
  bool chol_success = arma::inv_sympd(V_inv, V);
  if (!chol_success) {
    // Add small ridge
    V.diag() += 1e-6;
    arma::inv_sympd(V_inv, V);
  }
  
  // Fixed effects (BLUE): b = (X'V^{-1}X)^{-1} X'V^{-1}y
  arma::mat XtVinv = X.t() * V_inv;
  arma::mat XtVinvX = XtVinv * X;
  arma::vec XtVinvy = XtVinv * y;
  
  arma::vec b;
  arma::solve(b, XtVinvX, XtVinvy);
  
  // Compute P*y for BLUP
  arma::mat P = V_inv - V_inv * X * arma::solve(XtVinvX, XtVinv);
  arma::vec Py = P * y;
  
  // Random effects (BLUP): u = sigma2_g * K * Z' * P * y
  arma::vec u = sigma2_g * K * Z.t() * Py;
  
  // Fitted values
  arma::vec fitted = X * b + Z * u;
  
  // Residuals
  arma::vec residuals = y - fitted;
  
  return Rcpp::List::create(
    Rcpp::Named("b") = b,
    Rcpp::Named("u") = u,
    Rcpp::Named("fitted") = fitted,
    Rcpp::Named("residuals") = residuals,
    Rcpp::Named("Py") = Py
  );
}

//' Check if matrix is positive definite
//' 
//' @param M Square matrix
//' @param tol Tolerance for smallest eigenvalue
//' @return TRUE if positive definite
//' @export
// [[Rcpp::export]]
bool is_positive_definite_cpp(const arma::mat& M, double tol = 1e-8) {
  if (M.n_rows != M.n_cols) {
    return false;
  }
  
  // Ensure symmetry
  arma::mat M_sym = (M + M.t()) / 2.0;
  arma::vec eigval = arma::eig_sym(M_sym);
  return arma::min(eigval) > tol;
}

//' Make matrix positive definite by adding ridge
//' 
//' @param M Square matrix
//' @param tol Tolerance for smallest eigenvalue
//' @return Positive definite matrix
//' @export
// [[Rcpp::export]]
arma::mat make_pd_cpp(const arma::mat& M, double tol = 1e-6) {
  // Ensure input is symmetric
  arma::mat M_sym = (M + M.t()) / 2.0;
  
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, M_sym);
  
  // Floor eigenvalues
  for (size_t i = 0; i < eigval.n_elem; i++) {
    if (eigval(i) < tol) {
      eigval(i) = tol;
    }
  }
  
  // Reconstruct
  arma::mat M_pd = eigvec * arma::diagmat(eigval) * eigvec.t();
  
  // Ensure symmetry
  M_pd = (M_pd + M_pd.t()) / 2.0;
  
  return M_pd;
}

//' Fast PLINK BED File Decoding
//' 
//' Decodes PLINK binary BED format genotype data to integer matrix.
//' Uses standard PLINK BED encoding:
//' - 00 (binary) = 0 (homozygous reference)
//' - 01 (binary) = NA (missing)
//' - 10 (binary) = 1 (heterozygous)
//' - 11 (binary) = 2 (homozygous alternate)
//' 
//' @param bed_data Raw vector containing BED file bytes (after 3-byte header)
//' @param n_ind Number of individuals
//' @param n_snp Number of SNPs
//' @return Integer matrix (n_ind x n_snp) with genotypes coded as 0, 1, 2, or NA
//' @export
// [[Rcpp::export]]
IntegerMatrix decode_bed_cpp(RawVector bed_data, int n_ind, int n_snp) {
  // PLINK BED format stores SNP-major: each SNP's genotypes packed in bytes
  // Each byte contains up to 4 genotypes (2 bits each)
  
  IntegerMatrix geno(n_ind, n_snp);
  
  // Lookup table for 2-bit genotype decoding
  // PLINK encoding: 00=hom1, 01=missing, 10=het, 11=hom2
  int lookup[4] = {0, NA_INTEGER, 1, 2};
  
  // Calculate bytes per SNP (need ceiling division)
  int bytes_per_snp = (n_ind + 3) / 4;
  
  // Iterate over SNPs
  for (int snp = 0; snp < n_snp; snp++) {
    int byte_start = snp * bytes_per_snp;
    
    // Iterate over individuals
    for (int ind = 0; ind < n_ind; ind++) {
      // Calculate which byte and which position within byte
      int byte_idx = byte_start + (ind / 4);
      int bit_pos = (ind % 4) * 2;  // Each genotype uses 2 bits
      
      // Check bounds
      if (byte_idx >= bed_data.size()) {
        Rcpp::stop("BED data index out of bounds");
      }
      
      // Extract 2-bit genotype
      unsigned char byte_val = bed_data[byte_idx];
      int geno_code = (byte_val >> bit_pos) & 3;  // Extract 2 bits
      
      // Decode using lookup table
      geno(ind, snp) = lookup[geno_code];
    }
  }
  
  return geno;
}

//' Henderson's MME Solver in C++ (Single Random Effect)
//' 
//' Solves Henderson's Mixed Model Equations in (p+q)-space for maximum efficiency.
//' Much faster than observation-space methods when n >> q.
//' 
//' @param y Response vector (n x 1)
//' @param X Fixed effects design matrix (n x p)
//' @param Z Random effects design matrix (n x q)
//' @param K_inv Inverse of relationship matrix (q x q)
//' @param sigma2_g Genetic variance component
//' @param sigma2_e Residual variance component
//' @return List with fixed effects (b), random effects (u), fitted values, residuals
//' @export
// [[Rcpp::export]]
Rcpp::List solve_henderson_mme_cpp(const arma::vec& y, 
                                    const arma::mat& X,
                                    const arma::mat& Z, 
                                    const arma::mat& K_inv,
                                    double sigma2_g, 
                                    double sigma2_e) {
  
  int n = y.n_elem;
  int p = X.n_cols;
  int q = Z.n_cols;
  int dim = p + q;
  
  // Pre-compute cross products (efficient for sparse Z)
  arma::mat XtX = X.t() * X;
  arma::vec Xty = X.t() * y;
  arma::mat ZtZ = Z.t() * Z;
  arma::mat ZtX = Z.t() * X;
  arma::vec Zty = Z.t() * y;
  
  // Build Henderson's coefficient matrix C (Convention 2: R⁻¹-scaled)
  // C = [X'X/σ²_e        X'Z/σ²_e            ]
  //     [Z'X/σ²_e        Z'Z/σ²_e + K⁻¹/σ²_g ]
  // NOTE: Previously used lambda*K_inv where lambda=σ²_e/σ²_g (Convention 1 ratio),
  // but data blocks are already divided by σ²_e (Convention 2), creating a mixed
  // convention that makes the G⁻¹ penalty σ²_e times too strong.
  // FIX: Use K_inv/sigma2_g for consistent Convention 2.
  arma::mat C(dim, dim, arma::fill::zeros);
  
  // Top-left block: X'X / sigma2_e
  C.submat(0, 0, p-1, p-1) = XtX / sigma2_e;
  
  // Top-right and bottom-left blocks: X'Z / sigma2_e and Z'X / sigma2_e
  C.submat(0, p, p-1, dim-1) = ZtX.t() / sigma2_e;
  C.submat(p, 0, dim-1, p-1) = ZtX / sigma2_e;
  
  // Bottom-right block: Z'Z / sigma2_e + K_inv / sigma2_g (Convention 2)
  C.submat(p, p, dim-1, dim-1) = ZtZ / sigma2_e + K_inv / sigma2_g;
  
  // Build RHS vector
  arma::vec rhs(dim);
  rhs.subvec(0, p-1) = Xty / sigma2_e;
  rhs.subvec(p, dim-1) = Zty / sigma2_e;
  
  // Solve C [b; u] = rhs using Cholesky (faster and numerically stable)
  arma::vec sol;
  bool success = arma::solve(sol, C, rhs, arma::solve_opts::fast);
  
  if (!success) {
    // Fallback: add ridge
    C.diag() += DEFAULT_RIDGE;
    arma::solve(sol, C, rhs);
  }
  
  // Extract solutions
  arma::vec b = sol.subvec(0, p-1);
  arma::vec u = sol.subvec(p, dim-1);
  
  // Compute fitted values and residuals
  arma::vec fitted = X * b + Z * u;
  arma::vec residuals = y - fitted;
  
  return Rcpp::List::create(
    Rcpp::Named("b") = b,
    Rcpp::Named("u") = u,
    Rcpp::Named("fitted") = fitted,
    Rcpp::Named("residuals") = residuals,
    Rcpp::Named("C_inv") = arma::inv_sympd(C),  // For inference
    Rcpp::Named("C") = C,        // For exact REML LL computation
    Rcpp::Named("rhs") = rhs,    // For y'Py via Henderson identity
    Rcpp::Named("solution") = sol // For y'Py via Henderson identity
  );
}

//' Henderson's MME Solver with EM-REML (C++ Version)
//' 
//' Complete Henderson's MME solver with EM-REML variance component estimation.
//' Supports single random effect for now (can be extended to multiple).
//' 
//' @param y Response vector (n x 1)
//' @param X Fixed effects design matrix (n x p)
//' @param Z Random effects design matrix (n x q)
//' @param K Relationship matrix (q x q)
//' @param sigma2_g_init Initial genetic variance
//' @param sigma2_e_init Initial residual variance
//' @param max_iter Maximum iterations
//' @param tol Convergence tolerance
//' @param verbose Print progress
//' @return List with model components including variance components
//' @export
// [[Rcpp::export]]
Rcpp::List fit_henderson_mme_cpp(const arma::vec& y, 
                                  const arma::mat& X,
                                  const arma::mat& Z, 
                                  const arma::mat& K,
                                  double sigma2_g_init = -1.0,
                                  double sigma2_e_init = -1.0,
                                  int max_iter = 50,
                                  double tol = 1e-4,
                                  bool verbose = false) {
  
  int n = y.n_elem;
  int p = X.n_cols;
  int q = Z.n_cols;
  
  // Initialize variance components if not provided
  double var_y = arma::var(y);
  double sigma2_g = (sigma2_g_init > 0) ? sigma2_g_init : 0.5 * var_y;
  double sigma2_e = (sigma2_e_init > 0) ? sigma2_e_init : 0.5 * var_y;
  
  // Invert K once (q x q matrix)
  arma::mat K_inv;
  bool inv_success = arma::inv_sympd(K_inv, K);
  if (!inv_success) {
    // Add ridge if singular
    arma::mat K_ridge = K;
    K_ridge.diag() += DEFAULT_RIDGE;
    arma::inv_sympd(K_inv, K_ridge);
  }
  
  double llik_old = -arma::datum::inf;
  bool converged = false;
  int iter;
  
  arma::vec b, u, fitted, residuals;
  
  for (iter = 1; iter <= max_iter; iter++) {
    
    // Solve Henderson's MME with current variance components
    Rcpp::List sol = solve_henderson_mme_cpp(y, X, Z, K_inv, sigma2_g, sigma2_e);
    b = Rcpp::as<arma::vec>(sol["b"]);
    u = Rcpp::as<arma::vec>(sol["u"]);
    residuals = Rcpp::as<arma::vec>(sol["residuals"]);
    
    // EM updates for variance components with trace correction
    // BUG FIX (v2.3.1): Full EM-REML for sigma2_e with trace correction.
    // Previously: sigma2_e_new = rss / (n - p) (simplified, biased for multi-effect)
    double rss = arma::accu(arma::square(residuals));
    
    // Get C_inv from the solver (already computed for inference)
    arma::mat C_inv = Rcpp::as<arma::mat>(sol["C_inv"]);
    
    // Full EM-REML: σ²_g = (u'K⁻¹u + tr(C⁻¹_uu K⁻¹)) / q  (Convention 2)
    double uKu = arma::as_scalar(u.t() * K_inv * u);
    
    // Extract C⁻¹_uu block (rows/cols p to p+q-1)
    arma::mat C_inv_uu = C_inv.submat(p, p, p + q - 1, p + q - 1);
    
    // Trace correction for sigma2_g
    double trace_term = arma::trace(C_inv_uu * K_inv);
    double sigma2_g_new = std::max((uKu + trace_term) / q, MIN_VARIANCE);
    
    // Full EM-REML for sigma2_e with trace correction
    // tr(C⁻¹ · C_data) = (p+q) - tr(C⁻¹_uu · K⁻¹) / σ²_g
    double trace_Cinv_Ginv = trace_term / sigma2_g;
    double trace_Cinv_Cdata = (p + q) - trace_Cinv_Ginv;
    double sigma2_e_new = std::max((rss + sigma2_e * trace_Cinv_Cdata) / (n - p),
                                   MIN_VARIANCE);
    
    // BUG FIX (v2.3.1): Compute exact REML LL for convergence monitoring.
    // Previously used simplified -0.5*[n*log(σ²_e) + RSS/σ²_e].
    arma::mat C_mat = Rcpp::as<arma::mat>(sol["C"]);
    double log_det_C_val;
    double log_det_C_sign;
    arma::log_det(log_det_C_val, log_det_C_sign, C_mat);
    
    double log_det_K_val;
    double log_det_K_sign;
    arma::log_det(log_det_K_val, log_det_K_sign, K);
    
    double log_det_G = q * std::log(sigma2_g) + log_det_K_val;
    
    arma::vec rhs_vec = Rcpp::as<arma::vec>(sol["rhs"]);
    arma::vec sol_vec = Rcpp::as<arma::vec>(sol["solution"]);
    double yPy = arma::accu(arma::square(y)) / sigma2_e - arma::dot(rhs_vec, sol_vec);
    
    double llik = -0.5 * ((n - p) * std::log(sigma2_e) + log_det_G + log_det_C_val + yPy);
    
    // Check convergence
    double llik_diff = std::abs(llik - llik_old);
    double theta_diff = std::max(
      std::abs(sigma2_g_new - sigma2_g) / std::max(sigma2_g, 1e-8),
      std::abs(sigma2_e_new - sigma2_e) / std::max(sigma2_e, 1e-8)
    );
    
    if (llik_diff < tol && theta_diff < tol && iter > 1) {
      converged = true;
      if (verbose) {
        Rcpp::Rcout << "Converged at iteration " << iter << std::endl;
      }
      break;
    }
    
    // Update variance components
    sigma2_g = sigma2_g_new;
    sigma2_e = sigma2_e_new;
    llik_old = llik;
    
    if (verbose && iter % 5 == 0) {
      Rcpp::Rcout << "Iteration " << iter << ": LogLik = " << llik 
                  << ", sigma2_g = " << sigma2_g 
                  << ", sigma2_e = " << sigma2_e << std::endl;
    }
  }
  
  // Final solution with converged variance components
  Rcpp::List final_sol = solve_henderson_mme_cpp(y, X, Z, K_inv, sigma2_g, sigma2_e);
  b = Rcpp::as<arma::vec>(final_sol["b"]);
  u = Rcpp::as<arma::vec>(final_sol["u"]);
  fitted = Rcpp::as<arma::vec>(final_sol["fitted"]);
  residuals = Rcpp::as<arma::vec>(final_sol["residuals"]);
  
  // Compute final log-likelihood
  double rss = arma::accu(arma::square(residuals));
  double llik = -0.5 * (n * std::log(sigma2_e) + rss / sigma2_e);
  
  // Create theta vector
  arma::vec theta(2);
  theta(0) = sigma2_g;
  theta(1) = sigma2_e;
  
  return Rcpp::List::create(
    Rcpp::Named("b") = b,
    Rcpp::Named("u") = u,
    Rcpp::Named("fitted") = fitted,
    Rcpp::Named("residuals") = residuals,
    Rcpp::Named("sigma2_g") = sigma2_g,
    Rcpp::Named("sigma2_e") = sigma2_e,
    Rcpp::Named("theta") = theta,
    Rcpp::Named("llik") = llik,
    Rcpp::Named("convergence") = converged,
    Rcpp::Named("nIter") = iter,
    Rcpp::Named("C_inv") = final_sol["C_inv"]
  );
}
