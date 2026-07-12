#' Simulate Genotype Data
#'
#' Generates a simulated genotype matrix with specified dimensions and MAF range.
#'
#' @param n_ind Number of individuals
#' @param n_snp Number of SNPs
#' @param maf_range MAF range (default c(0.05, 0.5))
#' @return GBLUPData object with simulated genotypes
#' @export
#' @importFrom stats rbinom runif
#' @examples
#' \dontrun{
#' gdata <- simulate_genotypes(n_ind = 100, n_snp = 1000)
#' show(gdata)
#' }
simulate_genotypes <- function(n_ind, n_snp, maf_range = c(0.05, 0.5)) {
  
  # Validate inputs
  stopifnot("n_ind must be positive" = n_ind > 0)
  stopifnot("n_snp must be positive" = n_snp > 0)
  stopifnot("maf_range must have 2 elements" = length(maf_range) == 2)
  stopifnot("maf_range[1] must be less than maf_range[2]" = maf_range[1] < maf_range[2])
  stopifnot("maf_range must be between 0 and 0.5" = maf_range[1] >= 0 && maf_range[2] <= 0.5)
  
  # Generate MAF
  maf <- runif(n_snp, maf_range[1], maf_range[2])
  
  # Generate genotypes (0, 1, 2)
  G_raw <- matrix(
    rbinom(n_ind * n_snp, 2, rep(maf, each = n_ind)),
    nrow = n_ind, ncol = n_snp
  )
  
  # Convert to -1, 0, 1
  G_add <- G_raw - 1
  
  # Set names
  rownames(G_add) <- paste0("ID", seq_len(n_ind))
  colnames(G_add) <- paste0("SNP", seq_len(n_snp))
  
  # Create empty phenotype data frame
  pheno <- data.frame(
    GID = character(0),
    ENV = character(0),
    trait = numeric(0),
    stringsAsFactors = FALSE
  )
  
  # Return GBLUPData object
  new("GBLUPData",
      genotypes = G_add,
      phenotypes = pheno,
      maf = maf)
}

#' Simulate Phenotypes
#'
#' Simulates phenotypic values based on genotype data with specified genetic architecture.
#'
#' @param gblup_data GBLUPData object
#' @param n_env Number of environments
#' @param h2 Heritability (0-1)
#' @param n_qtl Number of causal QTL
#' @param include_dominance Include dominance effects
#' @param include_epistasis Include epistatic effects (AA)
#' @param include_gxe Include GxE effects
#' @param mu Overall mean
#' @param prop_dom Proportion of genetic variance due to dominance (default 0.2)
#' @param prop_epi Proportion of genetic variance due to epistasis (default 0.15)
#' @param prop_gxe Proportion of genetic variance due to GxE (default 0.1)
#' @param prop_env Proportion of environmental variance due to environment main effect (default 0.2)
#' @return GBLUPData object with simulated phenotypes
#' @export
#' @importFrom stats rnorm
#' @examples
#' \dontrun{
#' gdata <- simulate_genotypes(n_ind = 100, n_snp = 500)
#' gdata <- simulate_phenotypes(gdata, n_env = 3, h2 = 0.6, n_qtl = 20)
#' head(gdata@phenotypes)
#' }
simulate_phenotypes <- function(gblup_data, 
                                n_env = 3,
                                h2 = 0.6,
                                n_qtl = 50,
                                include_dominance = FALSE,
                                include_epistasis = FALSE,
                                include_gxe = TRUE,
                                mu = 100,
                                prop_dom = 0.2,
                                prop_epi = 0.15,
                                prop_gxe = 0.1,
                                prop_env = 0.2) {
  
  # Validate inputs
  stopifnot("gblup_data must be a GBLUPData object" = inherits(gblup_data, "GBLUPData"))
  stopifnot("n_env must be positive" = n_env > 0)
  stopifnot("h2 must be between 0 and 1" = h2 >= 0 && h2 <= 1)
  stopifnot("prop_dom must be between 0 and 1" = prop_dom >= 0 && prop_dom <= 1)
  stopifnot("prop_epi must be between 0 and 1" = prop_epi >= 0 && prop_epi <= 1)
  stopifnot("prop_gxe must be between 0 and 1" = prop_gxe >= 0 && prop_gxe <= 1)
  stopifnot("prop_env must be between 0 and 1" = prop_env >= 0 && prop_env <= 1)
  
  # Extract genotypes
  G <- gblup_data@genotypes
  n_ind <- nrow(G)
  n_snp <- ncol(G)
  
  stopifnot("n_qtl must not exceed n_snp" = n_qtl <= n_snp)
  
  # Select QTL
  qtl_idx <- sample(n_snp, n_qtl)
  
  # Simulate additive effects
  beta_add <- numeric(n_snp)
  beta_add[qtl_idx] <- rnorm(n_qtl, 0, 1)
  
  # Additive genetic values
  g_add <- as.numeric(G %*% beta_add)
  g_add <- scale(g_add) * sqrt(h2 * 100)  # Scale to h2
  
  # Initialize other effects
  g_dom <- 0
  g_epi <- 0
  
  # Dominance
  if (include_dominance) {
    G_dom <- (G == 0) * 1.0
    beta_dom <- numeric(n_snp)
    beta_dom[qtl_idx] <- rnorm(n_qtl, 0, 0.5)
    g_dom <- as.numeric(G_dom %*% beta_dom)
    g_dom <- scale(g_dom) * sqrt(h2 * 100 * prop_dom)
  }
  
  # Epistasis (AA)
  if (include_epistasis) {
    G_sq <- G * G
    beta_epi <- numeric(n_snp)
    beta_epi[qtl_idx] <- rnorm(n_qtl, 0, 0.3)
    g_epi <- as.numeric(G_sq %*% beta_epi)
    g_epi <- scale(g_epi) * sqrt(h2 * 100 * prop_epi)
  }
  
  # Total genetic value
  g_total <- g_add + g_dom + g_epi
  
  # Environmental variance
  var_g <- var(g_total)
  if (var_g == 0) var_g <- 1  # Avoid division by zero
  var_e <- var_g * (1 - h2) / h2
  
  # Generate phenotypes for each environment
  pheno_list <- lapply(seq_len(n_env), function(e) {
    
    # Environment main effect
    env_effect <- rnorm(1, 0, sqrt(var_e * prop_env))
    
    # GxE effect
    if (include_gxe) {
      beta_gxe <- rnorm(n_qtl, 0, 0.3)
      gxe_effect <- as.numeric(G[, qtl_idx] %*% beta_gxe)
      gxe_effect <- scale(gxe_effect) * sqrt(h2 * 100 * prop_gxe)
    } else {
      gxe_effect <- 0
    }
    
    # Residual error
    e_residual <- rnorm(n_ind, 0, sqrt(var_e * (1 - prop_env)))
    
    # Total phenotype
    y <- mu + g_total + gxe_effect + env_effect + e_residual
    
    data.frame(
      GID = rownames(G),
      ENV = paste0("E", e),
      trait = as.numeric(y),
      TBV = as.numeric(g_add),      # True ADDITIVE breeding value
      TGV = as.numeric(g_total),    # True TOTAL genetic value
      stringsAsFactors = FALSE
    )
  })
  
  # Combine all environments
  pheno_combined <- do.call(rbind, pheno_list)
  
  # Store simulation parameters for reference
  attr(pheno_combined, "sim_params") <- list(
    h2 = h2,
    n_qtl = n_qtl,
    qtl_idx = qtl_idx,
    include_dominance = include_dominance,
    include_epistasis = include_epistasis,
    include_gxe = include_gxe
  )
  
  # Update GBLUPData object
  gblup_data@phenotypes <- pheno_combined
  
  return(gblup_data)
}