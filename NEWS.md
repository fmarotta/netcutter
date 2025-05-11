# netcutter development version

# netcutter v0.3.0

* Introduce the `terms_of_interest` argument of nc_eval().
* Realize that the C++ speedup only worked for small-size matrices
* Devise algorithm that scales as (N + M) rather than (N * M), consequently drop Rcpp

# netcutter v0.2.1

* Rewrite the edge-swapping algorithm in C++ (it's now 20 times faster).
* Use combn() instead of expand.grid() in the internal function nc_define_modules().

As of today, the C++ version uses the old "Rounding" random sampler that R has used before version 3.6 ([see here](https://stackoverflow.com/questions/60119621/get-the-same-sample-of-integers-from-rcpp-as-base-r)).
This sampler is biased if the population size is large (some ppl say >1e7). Use v0.2.0 if you are concerned.

# netcutter v0.2.0

* Fix p-values calculation: provide `Pr(==k)` and `Pr(>=k)`.
* Add a `mc.cores` option to nc_eval()

The p-values calculation in previous releases was completely wrong: do not use releases older than v0.2.0!

# netcutter v0.1.1

* Fix bug occurring when setting R to 1 ([#1](https://github.com/fmarotta/netcutter/issues/1))
* Make nc_randomise() about two times faster
* Add a `verbose` option to nc_occ_probs()

# netcutter v0.1.0

* Initial release.
