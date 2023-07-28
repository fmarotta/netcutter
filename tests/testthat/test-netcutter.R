test_that("randomization works", {
  m <- matrix(F, nrow = 3, ncol = 9, dimnames = list(paste0("PMID", 1:3), paste0("gene", 1:9)))
  m[1, 1:3] <- m[2, c(1:2, 4:5)] <- m[3, c(1, 6:9)] <- T
  set.seed(1)
  m_random <- nc_randomize(m, 100)

  # the result should still be a logical matrix.
  expect_type(m_random, "logical")
  expect_equal(dim(m_random), dim(m))
  # the result should have the same number of edges.
  expect_equal(sum(m_random), sum(m))
  # the result should *probably* not be exactly identical to the starting matrix,
  # change the seed if this test fails.
  expect_false(identical(unname(m_random), unname(m)))
})

test_that("occurrence probabilities make sense", {
  m <- matrix(F, nrow = 3, ncol = 9, dimnames = list(paste0("PMID", 1:3), paste0("gene", 1:9)))
  m[1, 1:3] <- m[2, c(1:2, 4:5)] <- m[3, c(1, 6:9)] <- T

  # The complexity of mc.cores and n_batches shouldn't affect the end result
  set.seed(1, "L'Ecuyer-CMRG")
  occ_probs_ref <- nc_occ_probs_simple(m, R = 20, S = 50)
  set.seed(1, "L'Ecuyer-CMRG")
  occ_probs <- nc_occ_probs(m, R = 20, S = 50, mc.cores = 2)
  expect_equal(occ_probs, occ_probs_ref)
  # The complexity of mc.cores and n_batches shouldn't affect the end result
  set.seed(1, "L'Ecuyer-CMRG")
  occ_probs_ref <- nc_occ_probs_simple(m, R = 60, S = 50)
  set.seed(1, "L'Ecuyer-CMRG")
  occ_probs <- nc_occ_probs(m, R = 60, S = 50, mc.cores = 2, n_batches = 2)
  expect_equal(occ_probs, occ_probs_ref)
  set.seed(1, "L'Ecuyer-CMRG")
  occ_probs <- nc_occ_probs(m, R = 60, S = 50, mc.cores = 1, n_batches = 1)
  expect_equal(occ_probs, occ_probs_ref)
  # The probabilities should be between 0 and 1
  expect_true(all(occ_probs >= 0 & occ_probs <= 1))
  # The rows should add up to the number of edges of each container
  expect_identical(rowSums(occ_probs), rowSums(m))
  # The cols should add up to the number of edges of each term
  expect_identical(colSums(occ_probs), colSums(m))
})
