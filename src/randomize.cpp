#include <Rcpp.h>
#include <set>
using namespace Rcpp;
using namespace std;

// TODO: use a sparse matrix for m and mask.
// maybe a custom structure where all the elements in a row AND in a column are easily accessible.

// [[Rcpp::export]]
LogicalMatrix randomize(LogicalMatrix occ_matrix, unsigned int S) {
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
  NumericVector pool(n);
  int source, target;
  int c, r;
  int i, j, k, l;
  int mask_size;
  int m_nrow = m.nrow();
  int m_size = m.size();
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
    r = source % m_nrow;
    c = source / m_nrow;
    // We have to sample a source and a target such that they are the vertices
    // of a rectangle whose other two vertices are FALSE. In other words, after
    // we choose the source, the target cannot be in a row where the source's
    // column has true, and it cannot be in a column where the source's row has
    // true.
    // Iterating in this weird way is almost two times faster than iterating
    // by rows and columns with m(i, j)
    // Rcout << "source: " << r+1 << " " << c + 1 << endl;
    for (k = c * m_nrow, i = 0; i < m_nrow; ++k, ++i) {
      if (m[k]) {
        for (l = i; l < m_size; l += m_nrow) {
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
    for (k = r, j = 0; k < m_size; k += m_nrow, ++j) {
      if (m[k]) {
        i = j * m_nrow;
        for (l = i; l < i + m_nrow; ++l) {
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
    // Rcout << "target: " << target % m_nrow + 1 << " " << target / m_nrow + 1 << endl;
    // Swap
    m[(target / m_nrow) * m_nrow + r] = m[c * m_nrow + target % m_nrow] = true;
    m[source] = m[target] = false;
    s += 1;
  }
  return m;
}


/*** R
set.seed(2, "L'Ecu", sample.kind = "Rounding")
mr <- nc_randomize_R(m, 100)
set.seed(2, "L'Ecu", sample.kind = "Rounding")
mc <- randomize(m, 100)
identical(unname(mc), unname(mr))
microbenchmark(nc_randomize_R(m, 1000), randomize(m, 1000))
*/
