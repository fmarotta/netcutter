test_that("randomization works", {
  m <- matrix(FALSE, nrow = 3, ncol = 9, dimnames = list(paste0("PMID", 1:3), paste0("gene", 1:9)))
  m[1, 1:3] <- m[2, c(1:2, 4:5)] <- m[3, c(1, 6:9)] <- TRUE
  set.seed(1)
  m_random <- nc_randomize_fast(m, 100)

  # the result should still be a logical matrix.
  expect_type(m_random, "logical")
  expect_equal(dim(m_random), dim(m))
  # the result should have the same number of edges.
  expect_equal(sum(m_random), sum(m))
  # the result should *probably* not be exactly identical to the starting matrix,
  # change the seed if this test fails.
  expect_false(identical(unname(m_random), unname(m)))

  # the result should not depend on the implementation (within statistical fluctuations)
  R <- 500
  S <- 100
  set.seed(1, "L'Ecuyer-CMRG")
  seeds <- generate_seeds(R)
  swaps <- lapply(seeds, function(r) {
    .Random.seed <<- r
    nc_randomize_fast(m, S)
  })
  swaps_R <- lapply(seeds, function(r) {
    .Random.seed <<- r
    nc_randomize_R(m, S)
  })
  swaps_mean <- Reduce(`+`, swaps) / R
  swaps_mean_R <- Reduce(`+`, swaps_R) / R
  # expect_equal(unname(swaps_mean), swaps_mean_R)
  swaps_simple <- lapply(seeds, function(r) {
    .Random.seed <<- r
    nc_randomize_simple(m, S)
  })
  swaps_mean_simple <- Reduce(`+`, swaps_simple) / R
  # Swaps are 0-1 variables so their are equal to their square
  swaps_var_simple <- swaps_mean_simple - swaps_mean_simple^2
  expect_true(all(abs(swaps_mean - swaps_mean_simple) <= sqrt(swaps_var_simple)))
})

test_that("occurrence probabilities make sense", {
  m <- matrix(FALSE, nrow = 3, ncol = 9, dimnames = list(paste0("PMID", 1:3), paste0("gene", 1:9)))
  m[1, 1:3] <- m[2, c(1:2, 4:5)] <- m[3, c(1, 6:9)] <- TRUE

  # The complexity of n_batches shouldn't affect the end result
  set.seed(1, "L'Ecuyer-CMRG")
  occ_probs_ref <- nc_occ_probs_simple(m, R = 60, S = 50)
  set.seed(1, "L'Ecuyer-CMRG")
  occ_probs <- nc_occ_probs(m, R = 60, S = 50, n_batches = 3)
  expect_equal(occ_probs, occ_probs_ref)
  set.seed(1, "L'Ecuyer-CMRG")
  occ_probs <- nc_occ_probs(m, R = 60, S = 50, n_batches = 1)
  expect_equal(occ_probs, occ_probs_ref)
  # The probabilities should be between 0 and 1
  expect_true(all(occ_probs >= 0 & occ_probs <= 1))
  # The rows should add up to the number of edges of each container
  expect_identical(rowSums(occ_probs), rowSums(m))
  # The cols should add up to the number of edges of each term
  expect_identical(colSums(occ_probs), colSums(m))
  # Parameters should make sense
  expect_error(nc_occ_probs(m, R = 0))
  expect_error(nc_occ_probs(m, mc.cores = 0))
})

test_that("we can list all possible modules", {
  m <- matrix(FALSE, nrow = 3, ncol = 9, dimnames = list(paste0("PMID", 1:3), paste0("gene", 1:9)))
  m[1, 1:3] <- m[2, c(1:2, 4:5)] <- m[3, c(1, 6:9)] <- TRUE

  # With a size of 2, we expect choose(ncol(m), 2) possible modules.
  modules_two <- nc_define_modules(m, terms_of_interest = NULL, module_size = 2, min_occurrences = 0)
  expect_equal(length(modules_two), choose(ncol(m), 2))
  # With a size of 3, we expect choose(ncol(m), 3) possible modules.
  modules_three <- nc_define_modules(m, terms_of_interest = NULL, module_size = 3, min_occurrences = 0)
  expect_equal(length(modules_three), choose(ncol(m), 3))
  # All the terms should be represented
  expect_true(all(colnames(m) %in% unique(unlist(modules_two))))
  expect_true(all(unique(unlist(modules_two)) %in% colnames(m)))
  # Check that `min_occurrences` works
  modules <- nc_define_modules(m, terms_of_interest = NULL, module_size = 2, min_occurrences = 2)
  expect_equal(length(modules), choose(sum(colSums(m) >= 2), 2))
  modules <- nc_define_modules(m, terms_of_interest = NULL, module_size = 3, min_occurrences = 1)
  expect_equal(length(modules), choose(sum(colSums(m) >= 1), 3))
  # Check that `terms_of_interest` are included
  toi <- c("gene1", "gene9")
  modules <- nc_define_modules(m, terms_of_interest = toi, module_size = 2, min_occurrences = 0)
  expect_true(all(sapply(modules, function(x) any(toi %in% x))))
})

test_that("we can compute co-occurrence probabilities", {
  m <- matrix(FALSE, nrow = 3, ncol = 9, dimnames = list(paste0("PMID", 1:3), paste0("gene", 1:9)))
  m[1, 1:3] <- m[2, c(1:2, 4:5)] <- m[3, c(1, 6:9)] <- TRUE
  occ_probs <- nc_occ_probs(m, R = 20, S = 50)

  nc <- nc_eval(m, occ_probs, module_size = 2)
  # k should be between 0 and nrow(m)
  expect_true(all(nc$k >= 0 & nc$k < nrow(m)))
  # p-values should be between 0 and 1
  expect_true(all(nc$`Pr(==k)` >= 0 & nc$`Pr(==k)` <= 1))
  expect_true(all(nc$`Pr(>=k)` >= 0 & nc$`Pr(>=k)` <= 1))
  # If k < min_support, the p-values should be NA
  nc_minsupp <- nc_eval(m, occ_probs, module_size = 2, min_support = 1)
  expect_true(all(is.na(nc_minsupp[nc_minsupp$k < 1, ]$`Pr(==k)`)))
  expect_true(all(is.na(nc_minsupp[nc_minsupp$k < 1, ]$`Pr(>=k)`)))
})

test_that("generating random seeds works", {
  RNGkind("Mersenne")
  expect_error(generate_seeds(10))
  RNGkind("L'Ecu")
  expect_true(length(generate_seeds(10)) == 10)
  expect_true(length(generate_seeds(1)) == 1)
})

test_that("setting the seed works", {
  RNGkind("L'Ecu")
  test_rng <- function(r) {
    ans <- sample(5)
  }
  R <- 3
  set.seed(1)
  rlecuyer:::.lec.SetPackageSeed(sample(100, 6))
  rlecuyer:::.lec.CreateStream(seq_len(R))
  res <- parallel::mclapply(seq_len(R), mc.cores = 3, function(r) {
    oldkind <- rlecuyer:::.lec.CurrentStream(r)
    ans <- test_rng(r)
    rlecuyer:::.lec.CurrentStreamEnd(oldkind)
    ans
  })
  rlecuyer:::.lec.DeleteStream(seq_len(R))
})

test_that("parallelisation works (not on windows)", {
  skip_on_os("windows")

  m <- matrix(FALSE, nrow = 3, ncol = 9, dimnames = list(paste0("PMID", 1:3), paste0("gene", 1:9)))
  m[1, 1:3] <- m[2, c(1:2, 4:5)] <- m[3, c(1, 6:9)] <- TRUE

  # The complexity of mc.cores and n_batches shouldn't affect the end result
  set.seed(1, "L'Ecuyer-CMRG")
  occ_probs_ref <- nc_occ_probs_simple(m, R = 20, S = 50)
  set.seed(1, "L'Ecuyer-CMRG")
  occ_probs <- nc_occ_probs(m, R = 20, S = 50, mc.cores = 2)
  expect_equal(occ_probs, occ_probs_ref)
  set.seed(1, "L'Ecuyer-CMRG")
  occ_probs <- nc_occ_probs(m, R = 20, S = 50, mc.cores = 2, n_batches = 5)
  expect_equal(occ_probs, occ_probs_ref)

  # Using many cores shouldn't change the nc_eval() result
  nc_ref <- nc_eval(m, occ_probs_ref, mc.cores = 1)
  nc_twocores <- nc_eval(m, occ_probs_ref, mc.cores = 2)
  expect_equal(nc_ref, nc_twocores)
  nc_ref <- nc_eval(m, occ_probs_ref, min_support = 1, min_occurrences = 1, mc.cores = 1)
  nc_twocores <- nc_eval(m, occ_probs_ref, min_support = 1, min_occurrences = 1, mc.cores = 2)
  expect_equal(nc_ref, nc_twocores)
})
