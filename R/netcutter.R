#' Randomise the co-ocurrence matrix
#'
#' @param occ_matrix The original co-occurrence matrix.
#' @param S The number of successful edge swaps to perform.
#'
#' @return A randomised copy of the co-occurrence matrix.
#'
#' @examples
nc_randomize <- function(occ_matrix, S) {
  # Create a copy of the original matrix
  m <- matrix(occ_matrix, nrow(occ_matrix), ncol(occ_matrix))
  for (s in seq_len(S)) {
    es <- sapply(sample(which(m), 2), arrayInd, dim(m))
    while (m[es[1, 2], es[2, 1]] || m[es[1, 1], es[2, 2]]) {
      es <- sapply(sample(which(m), 2), arrayInd, dim(m))
    }
    m[es[1, 1], es[2, 1]] <- m[es[1, 2], es[2, 2]] <- F
    m[es[1, 1], es[2, 2]] <- m[es[1, 2], es[2, 1]] <- T
  }
  m
}

#' Compute the occurrence probabilities
#'
#' Use the EdgeSwapping method to find the probability of occurrence of each
#' term in each container under the null hypothesis.
#'
#' @param occ_matrix The original co-occurrence matrix
#' @param R The number of randomisations to perform
#' @param S The number of successful edge swaps for each randomisation
#' @param mc.cores Number of parallel computations with mclapply() (set to 1 for serial execution)
#' @param n_batches Split the computation into `n_batches` to avoid excessive memory usage

#' @return The occurrence probability matrix.
#' @export
#'
#' @examples
nc_occ_probs <- function(occ_matrix, R = 500, S = sum(occ_matrix) * 20,
                         mc.cores = getOption("mc.cores", 2L), n_batches = ceiling(R / 30)) {
  seeds <- generate_seeds(R)
  batch_size <- ceiling(R / n_batches)
  swaps <- matrix(0, nrow(occ_matrix), ncol(occ_matrix))
  for (batch in seq_len(n_batches)) {
    b <- seq((batch - 1) * batch_size + 1, min(batch * batch_size, R))
    swaps_batch <- parallel::mclapply(seeds[b], mc.cores = mc.cores, function(r) {
      .Random.seed <<- r
      nc_randomize(occ_matrix, S)
    })
    swaps <- swaps + Reduce(`+`, swaps_batch)
  }
  rownames(swaps) <- rownames(occ_matrix)
  colnames(swaps) <- colnames(occ_matrix)
  swaps / R
}

nc_occ_probs_simple <- function(occ_matrix, R, S) {
  seeds <- generate_seeds(R)
  swaps <- lapply(seeds, function(r) {
    .Random.seed <<- r
    nc_randomize(occ_matrix, S)
  })
  occ_probs <- Reduce(`+`, swaps) / R
  rownames(occ_probs) <- rownames(occ_matrix)
  colnames(occ_probs) <- colnames(occ_matrix)
  occ_probs
}

nc_define_modules <- function(occ_matrix, module_size, min_occurrences, min_support) {
  valid_terms <- which(colSums(occ_matrix) >= min_occurrences)
  message("The following terms are below the occurrence threshold: ",
          paste(which(colSums(occ_matrix) < min_occurrences), collapse = ", "))
  l <- rep(list(valid_terms), module_size)
  M <- expand.grid(l)
  if (module_size > 1) {
    for (j in seq(2, module_size)) {
      M <- eval(parse(text = paste0("M[M$Var", j - 1, " < ", "M$Var", j, ", ]")))
    }
  }
  support <- apply(M, 1, function(module) {
    sum(rowSums(occ_matrix[, module]) == module_size)
  })
  message("The following modules are below the support threshold: ",
          paste(apply(M[support < min_support, ], 1, paste, collapse = "-"), collapse = ","))
  M <- M[support >= min_support, ]
  if (!is.null(colnames(occ_matrix))) {
    M <- as.data.frame(lapply(M, function(x) colnames(occ_matrix)[x]))
  }
  rownames(M) <- NULL
  M
}


nc_eval <- function(occ_matrix, occ_probs, module_size = 2, min_occurrences = 0, min_support = 0) {
  M <- nc_define_modules(occ_matrix, module_size, min_occurrences, min_support)
  params <- apply(M, 1, function(module) {
    pp <- apply(occ_probs[, module], 1, prod)
    c(k, ppbinom(k, pp, method = "Characteristic"))
  })
  M$k <- params[1, ]
  M$pvals <- params[2, ]
  M
}

generate_seeds <- function(n) {
  seeds <- vector("list", n)
  seeds[[1]] <- .Random.seed
  for (i in 2:n)
    seeds[[i]] <- parallel::nextRNGStream(seeds[[i-1]])
  seeds
}





