test_that("add_ridge function works correctly", {
  # Create test matrix
  A <- matrix(c(1, 0.5, 0.5, 1), 2, 2)
  
  # Add ridge - use CompGEBLUP namespace explicitly if needed
  A_ridge <- CompGEBLUP::add_ridge(A, ridge = 0.01)
  
  # Tests
  expect_equal(dim(A_ridge), dim(A))
  expect_true(all(diag(A_ridge) > diag(A)))
  expect_equal(diag(A_ridge), diag(A) + 0.01)
  
  # Off-diagonal should be unchanged
  expect_equal(A_ridge[1, 2], A[1, 2])
})