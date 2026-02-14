#' Load PLINK BED File
#'
#' Reads PLINK binary BED format genotype data using fast C++ decoding.
#' Supports standard PLINK 1.9 BED format (SNP-major mode).
#'
#' @details
#' PLINK BED files store genotypes in binary format with 2 bits per genotype:
#' \itemize{
#'   \item 00 (binary) = 0 (homozygous reference)
#'   \item 01 (binary) = NA (missing)
#'   \item 10 (binary) = 1 (heterozygous)
#'   \item 11 (binary) = 2 (homozygous alternate)
#' }
#'
#' The BED file must be in SNP-major mode (standard for PLINK 1.9).
#' Requires corresponding .bim and .fam files to determine dimensions.
#'
#' @param bed_file Path to .bed file
#' @param bim_file Path to .bim file (if NULL, inferred from bed_file)
#' @param fam_file Path to .fam file (if NULL, inferred from bed_file)
#' @param coding Genotype coding: "012" (default) for 0/1/2, or "additive" for -1/0/1
#' @param impute_missing If TRUE, impute missing values with column means
#' @param verbose Print progress messages
#'
#' @return List with:
#'   \item{geno}{Integer or numeric matrix (n_ind x n_snp) of genotypes}
#'   \item{fam}{Data frame with individual information from .fam file}
#'   \item{bim}{Data frame with SNP information from .bim file}
#'
#' @export
#' @examples
#' \dontrun{
#' # Load PLINK files
#' plink_data <- load_plink_bed("mydata.bed")
#' 
#' # Extract genotype matrix coded as -1, 0, 1 for CompGEBLUP
#' M <- plink_data$geno - 1  # Convert 0/1/2 to -1/0/1
#' 
#' # Build relationship matrix
#' A <- build_A_matrix(M)
#' }
load_plink_bed <- function(bed_file, 
                           bim_file = NULL, 
                           fam_file = NULL,
                           coding = c("012", "additive"),
                           impute_missing = FALSE,
                           verbose = TRUE) {
  
  # Match arguments
  coding <- match.arg(coding)
  
  # Validate input
  if (!file.exists(bed_file)) {
    stop("BED file not found: ", bed_file)
  }
  
  # Infer .bim and .fam file paths if not provided
  if (is.null(bim_file)) {
    bim_file <- sub("\\.bed$", ".bim", bed_file)
  }
  if (is.null(fam_file)) {
    fam_file <- sub("\\.bed$", ".fam", bed_file)
  }
  
  if (!file.exists(bim_file)) {
    stop("BIM file not found: ", bim_file)
  }
  if (!file.exists(fam_file)) {
    stop("FAM file not found: ", fam_file)
  }
  
  # Read .fam file to get number of individuals
  if (verbose) message("Reading .fam file...")
  fam <- utils::read.table(fam_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(fam) <- c("FID", "IID", "PID", "MID", "Sex", "Phenotype")
  n_ind <- nrow(fam)
  
  # Read .bim file to get number of SNPs
  if (verbose) message("Reading .bim file...")
  bim <- utils::read.table(bim_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(bim) <- c("CHR", "SNP", "CM", "BP", "A1", "A2")
  n_snp <- nrow(bim)
  
  if (verbose) {
    message(sprintf("Loading genotypes: %d individuals x %d SNPs", n_ind, n_snp))
  }
  
  # Read BED file
  bed_size <- file.info(bed_file)$size
  expected_size <- 3 + ceiling(n_ind / 4) * n_snp
  
  if (bed_size != expected_size) {
    warning(sprintf(
      "BED file size (%d bytes) does not match expected size (%d bytes). File may be corrupted.",
      bed_size, expected_size
    ))
  }
  
  # Read raw bytes
  bed_con <- file(bed_file, "rb")
  on.exit(close(bed_con))
  
  # Check magic number (first 3 bytes: 0x6C, 0x1B, 0x01 for SNP-major)
  magic <- readBin(bed_con, "raw", n = 3)
  if (magic[1] != as.raw(0x6C) || magic[2] != as.raw(0x1B)) {
    stop("Invalid BED file magic number. File may be corrupted.")
  }
  if (magic[3] != as.raw(0x01)) {
    stop("BED file is in individual-major mode. Only SNP-major mode is supported.")
  }
  
  # Read genotype data (skip 3-byte header)
  bed_data <- readBin(bed_con, "raw", n = bed_size - 3)
  
  # Decode using C++ function
  if (verbose) message("Decoding genotypes...")
  geno <- decode_bed_cpp(bed_data, n_ind, n_snp)
  
  # Set row and column names
  rownames(geno) <- fam$IID
  colnames(geno) <- bim$SNP
  
  # Impute missing values if requested
  if (impute_missing) {
    if (verbose) message("Imputing missing values with column means...")
    for (j in 1:ncol(geno)) {
      missing_idx <- is.na(geno[, j])
      if (any(missing_idx)) {
        col_mean <- mean(geno[, j], na.rm = TRUE)
        geno[missing_idx, j] <- round(col_mean)
      }
    }
  }
  
  # Convert to additive coding if requested
  if (coding == "additive") {
    if (verbose) message("Converting to additive coding (-1/0/1)...")
    geno <- geno - 1
  }
  
  if (verbose) message("Done!")
  
  # Return list with genotypes and metadata
  list(
    geno = geno,
    fam = fam,
    bim = bim
  )
}
