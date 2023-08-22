#' Randomize the occurrence matrix
#'
#' Apply an edge-swapping algorithm.
#'
#' @param occ_matrix The original occurrence matrix.
#' @param S The number of successful edge swaps to perform.
#'
#' @return A randomized copy of the occurrence matrix.
nc_randomize <- function(occ_matrix, S) {
  stopifnot(S >= 0)
  .Call(`_netcutter_randomize`, occ_matrix, S)
}

#' Randomize the occurrence matrix
#'
#' Old implementation in pure R, kept for testing purposes and for
#' reproducibility of old results.
#'
#' @inheritParams nc_randomize
nc_randomize_R <- function(occ_matrix, S) {
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
#' @param verbose Print a status message when starting every new batch.
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
    nc_randomize(occ_matrix, S)
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
#' @inheritParams nc_eval
#'
#' @return A list of the valid modules.
nc_define_modules <- function(occ_matrix, terms_of_interest, module_size, min_occurrences) {
  valid_terms <- unname(which(colSums(occ_matrix) >= min_occurrences))
  if (!is.null(terms_of_interest)) {
    if (is.character(terms_of_interest) && all(terms_of_interest %in% colnames(occ_matrix))) {
      terms_of_interest <- match(terms_of_interest, colnames(occ_matrix))
    } else if (!(is.numeric(terms_of_interest) && all(terms_of_interest %in% 1:ncol(occ_matrix)))) {
      stop("`terms_of_interest` must be valid column names or indices for `occ_matrix`.")
    }
    other_terms <- setdiff(valid_terms, terms_of_interest)
    M <- combn(terms_of_interest, module_size, simplify = F)
    for (ms in seq_len(module_size - 1)) {
      toi_combn <- combn(terms_of_interest, ms, simplify = F)
      other_combn <- combn(other_terms, module_size - ms, simplify = F)
      M <- c(M, unlist(lapply(toi_combn, function(toi) lapply(other_combn, function(o) c(toi, o))), recursive = F))
    }
  } else {
    M <- utils::combn(valid_terms, module_size, simplify = F)
  }
  if (!is.null(colnames(occ_matrix))) {
    M <- lapply(M, function(x) colnames(occ_matrix)[x])
  }
  M
}

#' Compute co-occurrence probabilities
#'
#' The main NetCutter function. It generates p-values for all the co-occurring
#' modules.
#'
#' @param occ_matrix The original occurrence matrix.
#' @param occ_probs The matrix of occurrence probabilities, as computed by
#'   [nc_occ_probs()].
#' @param terms_of_interest Vector of column names or indices representing the
#'   terms that should be included in the analysis.
#' @param module_size The number of terms that should be tested for
#'   co-occurrence.
#' @param min_occurrences Minimum number of occurrences of each term.
#' @param min_support Minimum number of occurrences of each module.
#' @param mc.cores Number of parallel computations with mclapply() (set to 1 for
#'   serial execution)
#'
#' @details
#' If `terms_of_interest` is `NULL`, all the terms in `occ_matrix` are used. If
#' it is not null, only modules containing at least one of these terms will be
#' considered. `min_occurrences` and `min_support` are still used to further
#' restrict the list of terms that are considered.
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
#' # Now consider only modules involving gene1 or gene2.
#' nc_eval(m, occ_probs, module_size = 3, terms_of_interest = c("gene1", "gene2"))
#'
#' @importFrom PoissonBinomial ppbinom dpbinom
#' @export
nc_eval <- function(occ_matrix, occ_probs, terms_of_interest = NULL,
                    module_size = 2, min_occurrences = 0, min_support = 0,
                    mc.cores = 1) {
  stopifnot(module_size >= 1)
  M <- nc_define_modules(occ_matrix, terms_of_interest, module_size, min_occurrences)
  params <- parallel::mclapply(M, mc.cores = mc.cores, function(module) {
    k <- sum(rowSums(occ_matrix[, module]) == module_size)
    if (k < min_support) {
      return(c(k, NA, NA))
    }
    pp <- apply(occ_probs[, module], 1, prod)
    # pp <- eval(parse(text = paste(paste0("occ_probs[, module[", 1:module_size, "]]"), collapse = " * ")))
    # pp <- occ_probs[, module[1]]; for (m in 2:length(module)) pp <- pp * occ_probs[, module[2]]
    c(k,
      dpbinom(k, pp, method = "Characteristic"),
      ppbinom(k - 1, pp, method = "Characteristic", lower.tail = F))
  })
  res <- data.frame(
    do.call(rbind, M),
    do.call(rbind, params)
  )
  names(res) <- c(paste0("Term", 1:module_size), "k", "Pr(==k)", "Pr(>=k)")
  res
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
