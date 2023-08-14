#include <Rcpp.h>
#include <set>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
LogicalMatrix randomize(NumericMatrix occ_matrix, unsigned int S) {
  // NOTE: We need to use RNGkind(sample.kind = "Rounding") to get the same results. Rcpp only uses rounding, R can also use Rejection, which is better.
  LogicalMatrix m(clone(occ_matrix));
  LogicalMatrix mask(m.nrow(), m.ncol());
  LogicalMatrix::iterator it;
  int n = 0;
  for (it = m.begin(); it != m.end(); ++it) {
    if (*it) {
      ++n;
    }
  }
  vector<size_t> pool(n);
  size_t source, target;
  size_t c, r;
  int i, j, k, l;
  int mask_size;
  int nrow = m.nrow();
  unsigned int s = 0;
  while (s < S) {
    copy(m.begin(), m.end(), mask.begin());
    for (it = mask.begin(), l = 0; it != mask.end(); ++it) {
      if (*it) {
        pool[l++] = it - mask.begin();
      }
    }
    mask_size = n;
    source = pool[sample(mask_size, 1)[0] - 1];
    r = source % nrow;
    c = source / nrow;
    // Rcout << "source: " << r+1 << " " << c + 1 << endl;
    for (k = c * nrow, i = 0; i < nrow; ++k, ++i) {
      if (m[k]) {
        for (l = i; l < m.size(); l += nrow) {
          if (mask[l]) {
            mask[l] = false;
            --mask_size;
          }
        }
      }
    }
    // It's quite common to have terms in all container, so we have a shortcut.
    if (!mask_size) {
      continue;
    }
    for (k = r, j = 0; k < m.size(); k += nrow, ++j) {
      if (m[k]) {
        i = j * nrow;
        for (l = i; l < i + nrow; ++l) {
          if (mask[l]) {
            mask[l] = false;
            --mask_size;
          }
        }
      }
    }
    if (!mask_size) {
      continue;
    }
    for (it = mask.begin(), l = 0; it != mask.end(); ++it) {
      if (*it) {
        pool[l++] = it - mask.begin();
      }
    }
    target = pool[sample(mask_size, 1)[0] - 1];
    // Rcout << "target: " << target % nrow + 1 << " " << target / nrow + 1 << endl;
    // Swap
    m(r, target / nrow) = m(target % nrow, c) = true;
    m[source] = m[target] = false;
    s += 1;
  }
  return m;
}


/*** R
set.seed(2, "L'Ecu", sample.kind = "Rounding")
mr <- nc_randomize(m, 100)
set.seed(2, "L'Ecu", sample.kind = "Rounding")
mc <- randomize(m, 100)
identical(unname(mc), mr)
microbenchmark(nc_randomize(m, 1000), randomize(m, 1000))
*/
