#include <Rcpp.h>
#include <cstdint>
#include <x86intrin.h>

using namespace Rcpp;
using namespace std;

inline unsigned nthset(uint64_t x, unsigned n) {
    return _tzcnt_u64(_pdep_u64(1ULL << n, x));
}

// [[Rcpp::export]]
LogicalMatrix randomize_bitfiddling(const LogicalMatrix occ_matrix, int S) {
  const int lane_size = sizeof(uint64_t) * 8;
  size_t m_nrow = occ_matrix.nrow();
  size_t m_ncol = occ_matrix.ncol();
  size_t m_nlanes = (m_ncol + lane_size - 1) / lane_size;
  // cout << m_nlanes << endl;
  // cout << lane_size << endl;
  vector<vector<uint64_t>> rows(m_nrow, vector<uint64_t>(m_nlanes, 0));
  vector<uint64_t> unrolled_rows(m_nrow * m_nlanes, 0);
  vector<int> row_sums(m_nrow);
  size_t n_elements = 0;
  for (auto it = occ_matrix.begin(); it != occ_matrix.end(); ++it) {
    // cout << &(*it) << endl;
    // cout << *it << endl;
    int index = it - occ_matrix.begin();
    int r = index % m_nrow;
    int c = index / m_nrow;
    int l = m_nlanes - 1 - c / lane_size;
    int wlc = c % lane_size;
    rows[r][l] += static_cast<uint64_t>(*it) << wlc;
    row_sums[r] += *it;
    n_elements += *it;
  }
  /*
  cout << "rows" << endl;
  for (int r = 0; r < m_nrow; ++r) {
    cout << "row " << r << endl;
    cout << "row sum = " << row_sums[r] << endl;
    for (auto it = rows[r].begin(); it != rows[r].end(); ++it) {
      cout << &(*it) << endl;
      cout << hex << *it << dec << endl;
    }
  }
  */


  int s = 0;
  int source_col, source_row;
  int source_lane, source_wlc;
  uint64_t tempa, tempb, tempc;
  vector<uint64_t> candidates;
  uint64_t target;
  while (s++ < S) {
    // Sample the source (e.g. (alpha, A))
    // NOTE: if we want, we can reproduce R's code by doing everytihng by col rather than by row... OR, transpose the matrix in R so that it's row major...
    // cout << "Sampling" << endl;
    source_col = sample(n_elements, 1)[0] - 1; // 0-based
    source_row = 0; // 0-based
    while (source_col >= row_sums[source_row]) {
      source_col -= row_sums[source_row];
      ++source_row;
    }
    // cout << "Source row and col (zero-based)" << endl;
    // cout << source_row << endl;
    // cout << source_col << endl;

    // Get the col index of the source
    source_lane = m_nlanes - 1;
    while (source_col >= __builtin_popcountll(rows[source_row][source_lane])) {
      source_col -= __builtin_popcountll(rows[source_row][source_lane]);
      --source_lane;
    }
    source_wlc = nthset(rows[source_row][source_lane], source_col);
    /*
    cout << "The bit is in lane " << source_lane << endl;
    cout << "We look for the " << source_col << "th bit set" << endl;
    cout << source_lane << " " << source_wlc << endl;
    cout << endl;
    */

    /*
    // Find the candidates
    candidates.erase();
    for (int r = 0; r < m_nrow; ++r) {
      if (((rows[r][l] >> wlc) & 1ULL) == 1ULL) {
        continue;
      }
      for (int new_l = m_nlanes - 1; new_l >= 0; --new_l) {
        uint64_t tempa = ~rows[source_row][new_l] & rows[r][new_l];
        while (tempa != 0) {
          uint64_t t = tempa & -tempa;
          int new_wlc = __builtin_ctzll(tempa);
          candidates.push_back(vector<int>{r, new_l, new_wlc});
          tempa ^= t;
        }
      }
    }
    */
    // branch-free version
    candidates.clear();
    for (uint64_t candidate_row = 0; candidate_row < m_nrow; ++candidate_row) {
      tempa = ((rows[candidate_row][source_lane] >> source_wlc) & 1ULL) ^ 1ULL;
      // cout << "row " << candidate_row << " tempa " << tempa << endl;
      for (uint64_t candidate_lane = 0; candidate_lane < m_nlanes; ++candidate_lane) {
        tempb = (~rows[source_row][candidate_lane] & rows[candidate_row][candidate_lane]) * tempa;
        while (tempb != 0) {
          tempc = tempb & -tempb;
          uint64_t candidate_wlc = __builtin_ctzll(tempb);
          // candidates.push_back(vector<int>{candidate_row, candidate_lane, candidate_wlc});
          // uint64_t tempd = (candidate_row << 42) + (candidate_lane << 21) + candidate_wlc;
          candidates.push_back((candidate_row << 42) + (candidate_lane << 21) + candidate_wlc);
          tempb ^= tempc;
        }
      }
    }
    /*
    cout << "Candidates" << endl;
    for (int r = 0; r < candidates.size(); ++r) {
      cout << hex << candidates[r] << dec << endl;
      cout << "row " << ((candidates[r] >> 42) & 0x00000000001fffff) << endl;
      cout << "lane " << ((candidates[r] >> 21) & 0x00000000001fffff) << endl;
      cout << "wlc " << (candidates[r] & 0x00000000001fffff) << endl;
    }
    */

    // Sample the candidate
    // cout << candidates.size() << endl;
    if (candidates.size() == 0) {
      continue;
    }
    uint64_t mask = 0x00000000001fffff;
    target = candidates[sample(candidates.size(), 1)[0] - 1];
    /*
    cout << "Target row " << ((target >> 42) & mask) << endl;
    cout << "Target lane " << ((target >> 21) & mask) << endl;
    cout << "Target wlc " << (target & mask) << endl;
    cout << "Source row " << source_row << endl;
    cout << "Source lane " << source_lane << endl;
    cout << "Source wlc " << source_wlc << endl;
    */

    // Swap
    rows[source_row][source_lane] ^= 1ULL << source_wlc;
    rows[source_row][(target >> 21) & mask] ^= 1ULL << (target & mask);
    rows[(target >> 42) & mask][source_lane] ^= 1ULL << source_wlc;
    rows[(target >> 42) & mask][(target >> 21) & mask] ^= 1ULL << (target & mask);
  }

  return occ_matrix;
}

/***R
m <- sapply(1:80, function(row) {
  nt <- sample(1:4, 1)
  row <- rep(F, 5)
  row[sample(5, nt)] <- T
  row
})
m <- matrix(0, 2, 6)
m[1, c(1, 4, 6)] <- 1
m[2, c(3, 4, 5)] <- 1
m <- matrix(3, 5)
randomize_bitfiddling(m, 1000)
*/

/*
// [[Rcpp::export]]
LogicalMatrix randomize2(LogicalMatrix occ_matrix, int S) {
  // NOTE: We need to use RNGkind(sample.kind = "Rounding") to get the same results. Rcpp only uses rounding, R can also use Rejection, which is better.
  LogicalMatrix m(clone(occ_matrix));
  LogicalMatrix mask(clone(m));
  int m_nrow = m.nrow();
  int m_ncol = m.ncol();
  int n = accumulate(m.begin(), m.end(), 0);
  int mask_size;
  NumericVector pool(n);
  vector<unordered_set<int>> rows(m_nrow);
  vector<unordered_set<int>> cols(m_ncol);
  int source, target;
  int c, r;
  int cc, rr;
  int k, l;
  for (auto it = mask.begin(); it != mask.end(); ++it) {
    if (*it) {
      r = (it - mask.begin()) % m_nrow;
      c = (it - mask.begin()) / m_nrow;
      rows[r].insert(c);
      cols[c].insert(r);
    }
  }
  int s = 0;
  while (s < S) {
    for (int j = 0, l = 0; j < m_ncol; ++j) {
      for (int i : cols[j]) {
        k = j * m_nrow + i;
        pool[l++] = k;
      }
    }
    mask_size = n;
    source = sample(pool, 1)[0];
    // source = pool[sample(mask_size, 1)[0] - 1];
    r = source % m_nrow;
    c = source / m_nrow;
    // We have to sample a source and a target such that they are the vertices
    // of a rectangle whose other two vertices are FALSE. In other words, after
    // we choose the source, the target cannot be in a row where the source's
    // column has true, and it cannot be in a column where the source's row has
    // true.
    // Rcout << "source: " << source << " (" << r+1 << ", " << c + 1 << ")" << endl;
    for (int j : rows[r]) {
      for (int i : cols[j]) {
        k = j * m_nrow + i;
        if (mask[k]) {
          mask[k] = false;
          --mask_size;
        }
      }
    }
    if (!mask_size) {
      for (int j = 0, l = 0; j < m_ncol; ++j) {
        for (int i : cols[j]) {
          k = j * m_nrow + i;
          if (!mask[k]) {
            mask[k] = true;
          }
        }
      }
      continue;
    }
    for (int i : cols[c]) {
      for (int j : rows[i]) {
        k = j * m_nrow + i;
        if (mask[k]) {
          mask[k] = false;
          --mask_size;
        }
      }
    }
    if (!mask_size) {
      for (int j = 0, l = 0; j < m_ncol; ++j) {
        for (int i : cols[j]) {
          k = j * m_nrow + i;
          if (!mask[k]) {
            mask[k] = true;
          }
        }
      }
      continue;
    }
    for (int j = 0, l = 0; j < m_ncol; ++j) {
      for (int i : cols[j]) {
        k = j * m_nrow + i;
        if (mask[k]) {
          pool[l++] = k;
        } else {
          mask[k] = true;
        }
      }
    }
    target = pool[sample(mask_size, 1)[0] - 1];
    rr = target % m_nrow;
    cc = target / m_nrow;
    // Rcout << "target: " << target << " (" << rr+1 << ", " << cc + 1 << ")" << endl;
    // Swap
    m[cc * m_nrow + r] = m[c * m_nrow + rr] = true;
    m[source] = m[target] = false;
    mask[cc * m_nrow + r] = mask[c * m_nrow + rr] = true;
    mask[source] = mask[target] = false;
    rows[r].erase(c);
    rows[r].insert(cc);
    rows[rr].erase(cc);
    rows[rr].insert(c);
    cols[c].erase(r);
    cols[c].insert(rr);
    cols[cc].erase(rr);
    cols[cc].insert(r);
    s += 1;
  }
  return m;
}
*/

/***R
m <- matrix(FALSE, 3, 9, dimnames = list(paste0("ID", 1:3), paste0("gene", 1:9)))
m[1, 1:3] <- m[2, c(1:2, 4:5)] <- m[3, c(1, 6:9)] <- TRUE
randomize2(m, 1)
*/


// [[Rcpp::export]]
LogicalMatrix randomize(LogicalMatrix occ_matrix, unsigned int S) {
  return occ_matrix;
}
/*
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
*/

/*** R
set.seed(20, "L'Ecu", sample.kind = "Rounding")
mr <- nc_randomize_R(m, 1)
set.seed(20, "L'Ecu", sample.kind = "Rounding")
mc2 <- randomize2(m, 1)
identical(unname(mc2), unname(mr))
set.seed(20, "L'Ecu", sample.kind = "Rounding")
mc <- randomize(m, 1)
identical(unname(mc), unname(mr))

set.seed(2, "L'Ecu", sample.kind = "Rounding")
mc <- randomize(m, 100)
identical(unname(mc), unname(mr))
identical(unname(mc2), unname(mc))

microbenchmark(nc_randomize_R(m, 1000), randomize2(m, 1000))

huge <- t(sapply(1:100, function(row) {
  nt <- sample(100, 1)
  row <- rep(F, 1000)
  row[sample(1000, nt)] <- T
  row
}))
microbenchmark(nc_randomize_R(huge, 100), randomize_bitfiddling(huge, 100), times = 10)

*/
