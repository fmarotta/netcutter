#' Randomize the occurrence matrix
#'
#' Apply an edge-swapping algorithm.
#'
#' @param occ_matrix The original occurrence matrix.
#' @param S The number of successful edge swaps to perform.
#'
#' @return A randomized copy of the occurrence matrix.
nc_randomize <- function(occ_matrix, S) {
  # Create a copy of the original matrix
  m <- matrix(occ_matrix, nrow(occ_matrix), ncol(occ_matrix))
  l <- length(m)
  d <- dim(m)
  s <- 1
  while (s <= S) {
    # Choose the first edge at random and the second among the "compatible ones"
    source <- arrayInd(sample((1:l)[m], 1), d)
    mask <- m
    mask[, which(m[source[1], ])] <- F
    mask[which(m[, source[2]]), ] <- F
    if (sum(mask) == 0) {
      # This source edge has no compatible targets
      next()
    }
    target <- arrayInd(sample((1:l)[mask], 1), d)
    m[source[1], source[2]] <- m[target[1], target[2]] <- F
    m[source[1], target[2]] <- m[target[1], source[2]] <- T
    s <- s + 1
  }
  m
}
nc_randomise <- nc_randomize

#' Randomize the occurrence matrix
#'
#' This is a simpler implementation used to check that the official
#' implementation ([nc_randomize()]) works well.
#'
#' @inheritParams nc_randomize
nc_randomize_simple <- function(occ_matrix, S) {
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
#' @param mc.cores Number of parallel computations with mclapply() (set to 1 for
#'   serial execution)
#' @param n_batches Split the computation into `n_batches` to avoid excessive
#'   memory usage
#'
#' @return The occurrence probability matrix.
#'
#' @examples
#' # Generate an occurrence matrix.
#' m <- matrix(FALSE, 3, 9, dimnames = list(paste0("ID", 1:3), paste0("gene", 1:9)))
#' m[1, 1:3] <- m[2, c(1:2, 4:5)] <- m[3, c(1, 6:9)] <- TRUE
#' # Set the seed using the "L'Ecuyer-CMRG" random number generator.
#' set.seed(1, "L'Ecuyer-CMRG")
#' # Compute the occurrence probabilities.
#' occ_probs <- nc_occ_probs(m, R = 20, S = 50)
#' # Using `n_batches=1` can speed up the computations at the cost of more RAM.
#' occ_probs <- nc_occ_probs(m, R = 20, n_batches = 1, mc.cores = 1)
#'
#' @export
nc_occ_probs <- function(occ_matrix, R = 500, S = sum(occ_matrix) * 10,
                         mc.cores = getOption("mc.cores", 1L), n_batches = ceiling(R / 30),
                         verbose = F) {
  batch_size <- ceiling(R / n_batches)
  if (!is.numeric(mc.cores) || mc.cores <= 0) {
    stop("mc.cores must be positive.")
  }
  if (batch_size < mc.cores) {
    stop("Because you are using ", mc.cores, " cores, the maximum number of ",
         "batches is ceiling(R / ", mc.cores, ") = ", ceiling(R / mc.cores))
  }
  seeds <- generate_seeds(R)
  swaps <- matrix(0, nrow(occ_matrix), ncol(occ_matrix))
  for (batch in seq_len(n_batches)) {
    if (verbose) {
      message("Starting batch ", batch)
    }
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

#' Compute the occurrence probabilities
#'
#' This is a simpler implementation used to check that the official
#' implementation ([nc_occ_probs()]) works well.
#'
#' @inheritParams nc_occ_probs
nc_occ_probs_simple <- function(occ_matrix, R, S) {
  seeds <- generate_seeds(R)
  swaps <- lapply(seeds, function(r) {
    .Random.seed <<- r
    nc_randomize_simple(occ_matrix, S)
  })
  occ_probs <- Reduce(`+`, swaps) / R
  rownames(occ_probs) <- rownames(occ_matrix)
  colnames(occ_probs) <- colnames(occ_matrix)
  occ_probs
}

#' Define co-occurrence modules
#'
#' Helper function to generate the list of co-occurrence terms grouped into
#' modules of a specified size.
#'
#' @param occ_matrix The original occurrence matrix.
#' @param module_size The number of terms that should be tested for
#'   co-occurrence.
#' @param min_occurrences Minimum number of occurrences of each term.
#'
#' @return A `data.frame` with one row for each valid module.
nc_define_modules <- function(occ_matrix, module_size, min_occurrences) {
  valid_terms <- which(colSums(occ_matrix) >= min_occurrences)
  l <- rep(list(valid_terms), module_size)
  M <- expand.grid(l)
  if (module_size > 1) {
    for (j in seq(2, module_size)) {
      M <- eval(parse(text = paste0("M[M$Var", j - 1, " < ", "M$Var", j, ", ]")))
    }
  }
  if (!is.null(colnames(occ_matrix))) {
    M <- as.data.frame(lapply(M, function(x) colnames(occ_matrix)[x]))
  }
  rownames(M) <- NULL
  M
}

#' Compute co-occurrence probabilities
#'
#' The main NetCutter function. It generates p-values for all the co-occurring
#' modules.
#'
#' @inheritParams nc_define_modules
#' @param occ_matrix The original occurrence matrix.
#' @param occ_probs The matrix of occurrence probabilities, as computed by
#'   [nc_occ_probs()].
#' @param min_support Minimum number of occurrences of each module.
#'
#' @return A `data.frame` with one row for each valid module, and corresponding
#'   number of co-occurrences and p-value.
#'
#' @examples
#' # Generate an occurrence matrix.
#' m <- matrix(FALSE, 3, 9, dimnames = list(paste0("ID", 1:3), paste0("gene", 1:9)))
#' m[1, 1:3] <- m[2, c(1:2, 4:5)] <- m[3, c(1, 6:9)] <- TRUE
#' # Set the seed using the "L'Ecuyer-CMRG" random number generator.
#' set.seed(1, "L'Ecuyer-CMRG")
#' # Compute the occurrence probabilities.
#' occ_probs <- nc_occ_probs(m, R = 20, S = 50)
#' # Evaluate the co-occurrences of pairs of terms and their statistical significance.
#' nc_eval(m, occ_probs, module_size = 2)
#' # Now evaluate triples; no need to recompute the occurrence probabilities.
#' nc_eval(m, occ_probs, module_size = 3)
#'
#' @importFrom PoissonBinomial ppbinom
#' @export
nc_eval <- function(occ_matrix, occ_probs, module_size = 2, min_occurrences = 0, min_support = 0) {
  M <- nc_define_modules(occ_matrix, module_size, min_occurrences)
  params <- apply(M, 1, function(module) {
    k <- sum(rowSums(occ_matrix[, module]) == module_size)
    if (k < min_support) {
      return(c(k, NA))
    }
    pp <- apply(occ_probs[, module], 1, prod)
    c(k, ppbinom(k, pp, method = "Characteristic"))
  })
  M$k <- params[1, ]
  M$pvals <- params[2, ]
  M
}

#' Generate seeds for streams of "L'Ecuyer-CMRG" random numbers
#'
#' @param n Number of streams needed.
#'
#' @returns A list of seeds.
generate_seeds <- function(n) {
  if (RNGkind()[1] != "L'Ecuyer-CMRG") {
    stop("Due to parallel computations and reproducibility, the netcutter",
         "package requires a \"L'Ecuyer-CMRG\" random number generator.",
         "Please run either `RNGKind(\"L'Ecuyer-CMRG\")` or",
         "`set.seed(<my seed>, 'L'Ecuyer-CMRG'.")
  }
  seeds <- vector("list", n)
  seeds[[1]] <- .Random.seed
  for (i in seq_len(n-1))
    seeds[[i+1]] <- parallel::nextRNGStream(seeds[[i]])
  seeds
}
