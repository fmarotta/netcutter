# netcutter

This is an R implementation of the NetCutter method described in [Identification and Analysis of Co-Occurrence Networks with NetCutter](https://doi.org/10.1371/journal.pone.0003178) by Heiko Müller and Francesco Mancuso.
Co-occurrence analysis is used to find groups of terms that occur together in more (or fewer) containers than you would expect by chance.
For instance, the terms could be words and the containers could be articles: two words would co-occur if they tend to occur together in many articles.
The approach used by NetCutter to evaluate the statistical significance of co-occurrences is to apply a randomisation procedure to the occurrence matrix to obtain the occurrence probabilities under the null distribution, then use the [Poisson-Binomial distribution](https://en.wikipedia.org/wiki/Poisson_binomial_distribution) to compute a p-value for the actual number of occurrences.
The randomisation procedure, based on edge-swapping, approximates the more rigorous approach, which requires generating all the edge permutations of the occurrence graph and is therefore computationally intractable.

## NOTE

I have made reasonable efforts to contact the original authors of the paper, but haven't heard from them so far.
They are also listed as authors of this R package because they came up with the original idea (all the credit goes to them).
I just liked the paper and, not finding any other implementation, decided to write an R package (all mistakes and bugs are due to me, not to the original authors).

## Installation

The pacakge can be installed from CRAN.

``` r
install.packages("netcutter")
```

You can install the development version of netcutter from [GitHub](https://github.com/fmarotta/netcutter) with:

``` r
# install.packages("devtools")
devtools::install_github("fmarotta/netcutter")
```

## Badgers

<!-- badges: start -->
[![R-CMD-check](https://github.com/fmarotta/netcutter/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fmarotta/netcutter/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/fmarotta/netcutter/graph/badge.svg)](https://app.codecov.io/gh/fmarotta/netcutter)
[![CRAN version](https://www.r-pkg.org/badges/version-ago/netcutter)](https://www.r-pkg.org/badges/version-ago/netcutter)
[![CRAN downloads](https://cranlogs.r-pkg.org/badges/netcutter)](https://cranlogs.r-pkg.org/badges/netcutter)

<!-- badges: end -->
