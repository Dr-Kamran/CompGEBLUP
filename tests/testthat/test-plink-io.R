test_that("decode_bed_cpp decodes PLINK BED format correctly", {
  # PLINK BED encoding:
  # 00 (binary) = 0 (homozygous reference)
  # 01 (binary) = NA (missing)
  # 10 (binary) = 1 (heterozygous)
  # 11 (binary) = 2 (homozygous alternate)
  
  # Test case 1: Simple 4 individuals x 2 SNPs
  # SNP 1: genotypes 0, NA, 1, 2 (binary: 00, 01, 10, 11 = 0xE4)
  # SNP 2: genotypes 2, 1, 0, NA (codes: 11, 10, 00, 01 → 01_00_10_11 = 0x4B)
  
  bed_data <- as.raw(c(0xE4, 0x4B))
  n_ind <- 4
  n_snp <- 2
  
  geno <- CompGEBLUP::decode_bed_cpp(bed_data, n_ind, n_snp)
  
  # Check dimensions
  expect_equal(nrow(geno), n_ind)
  expect_equal(ncol(geno), n_snp)
  
  # Check SNP 1 (0, NA, 1, 2)
  expect_equal(geno[1, 1], 0)
  expect_true(is.na(geno[2, 1]))
  expect_equal(geno[3, 1], 1)
  expect_equal(geno[4, 1], 2)
  
  # Check SNP 2 (2, 1, 0, NA)
  expect_equal(geno[1, 2], 2)
  expect_equal(geno[2, 2], 1)
  expect_equal(geno[3, 2], 0)
  expect_true(is.na(geno[4, 2]))
  
  
  # Test case 2: Edge case with 1 individual x 1 SNP
  # Genotype: 0 (binary: 00)
  bed_data_small <- as.raw(0x00)
  geno_small <- CompGEBLUP::decode_bed_cpp(bed_data_small, 1, 1)
  
  expect_equal(dim(geno_small), c(1, 1))
  expect_equal(geno_small[1, 1], 0)
  
  
  # Test case 3: Multiple individuals requiring multiple bytes
  # 5 individuals (2 bytes needed per SNP) x 1 SNP
  # Genotypes: 0, 1, 2, NA, 0
  # Byte 1: codes 00, 10, 11, 01 → 01_11_10_00 = 0x78
  # Byte 2: 00, XX, XX, XX = 0x00 (only first 2 bits used)
  
  bed_data_multi <- as.raw(c(0x78, 0x00))
  geno_multi <- CompGEBLUP::decode_bed_cpp(bed_data_multi, 5, 1)
  
  expect_equal(dim(geno_multi), c(5, 1))
  expect_equal(geno_multi[1, 1], 0)
  expect_equal(geno_multi[2, 1], 1)
  expect_equal(geno_multi[3, 1], 2)
  expect_true(is.na(geno_multi[4, 1]))
  expect_equal(geno_multi[5, 1], 0)
})


test_that("load_plink_bed validates inputs correctly", {
  skip_if_not(file.exists("test_data/sample.bed"), 
              "Test PLINK files not available")
  
  # Test missing file
  expect_error(
    CompGEBLUP::load_plink_bed("nonexistent.bed"),
    "BED file not found"
  )
  
  # If test files exist, test loading
  if (file.exists("test_data/sample.bed")) {
    result <- CompGEBLUP::load_plink_bed(
      "test_data/sample.bed",
      coding = "012",
      verbose = FALSE
    )
    
    expect_type(result, "list")
    expect_true("geno" %in% names(result))
    expect_true("fam" %in% names(result))
    expect_true("bim" %in% names(result))
    expect_true(is.matrix(result$geno))
    expect_true(is.data.frame(result$fam))
    expect_true(is.data.frame(result$bim))
  }
})
