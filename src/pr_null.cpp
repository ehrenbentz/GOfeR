#include <Rcpp.h>
#include <algorithm>
#include <vector>
#include <random>
#include <cstdint>
using namespace Rcpp;

// Random sets come from a 64-bit Mersenne Twister (std::mt19937_64, whose
// output sequence is fixed by the C++ standard) seeded per chunk, so results
// are identical on every platform and R's random number state is not used.
// A random set is the first n entries of a random permutation of 1..N
// (partial Fisher-Yates shuffle); unbiased integers by rejection sampling.
static inline int uniform_below(std::mt19937_64& rng, uint64_t m) {
    const uint64_t limit = UINT64_MAX - (UINT64_MAX % m);
    uint64_t x;
    do { x = rng(); } while (x >= limit);
    return (int) (x % m);
}

// The first max_n entries of B random permutations of 1..N (one row each).
// [[Rcpp::export(name = ".random_sets", rng = false)]]
IntegerMatrix random_sets(int B, int N, int max_n, double seed) {
    std::mt19937_64 rng((uint64_t) seed);
    std::vector<int> perm(N);
    for (int i = 0; i < N; i++) perm[i] = i + 1;
    std::vector<int> swaps(max_n);
    IntegerMatrix out(B, max_n);
    for (int b = 0; b < B; b++) {
        for (int n = 0; n < max_n; n++) {
            int k = n + uniform_below(rng, (uint64_t) (N - n));
            swaps[n] = k;
            std::swap(perm[n], perm[k]);
            out(b, n) = perm[n];
        }
        for (int n = max_n - 1; n >= 0; n--) std::swap(perm[n], perm[swaps[n]]);
    }
    return out;
}

// Average precision of B random gene sets for every category size at once.
// Ranks of a random set are inserted one at a time; two Fenwick trees hold
// the count and the sum of 1/rank of the inserted ranks. With
// S = sum_k k / r_(k) over the sorted set, inserting rank r adds (number of
// earlier ranks below r + 1) / r plus the sum of 1/r_j over earlier ranks
// above r. The average precision of the first n ranks is S / n.
// col_of[n - 1] is the output column (1-based) for set size n, 0 if size n
// is not needed.
// [[Rcpp::export(name = ".ap_null_chunk", rng = false)]]
NumericMatrix ap_null_chunk(int B, int N, IntegerVector col_of, int n_sizes, double seed) {
    int max_n = col_of.size();
    NumericMatrix out(B, n_sizes);
    std::mt19937_64 rng((uint64_t) seed);
    std::vector<int> perm(N), swaps(max_n);
    for (int i = 0; i < N; i++) perm[i] = i + 1;
    // Fenwick nodes carry the index of the random set that last wrote them;
    // a node written by an earlier set counts as zero, so the trees need no
    // reset between sets.
    std::vector<int> cnt(N + 1), stamp(N + 1, -1);
    std::vector<double> inv(N + 1);
    for (int b = 0; b < B; b++) {
        double S = 0.0, total_inv = 0.0;
        for (int n = 0; n < max_n; n++) {
            int k = n + uniform_below(rng, (uint64_t) (N - n));
            swaps[n] = k;
            std::swap(perm[n], perm[k]);
            int r = perm[n];
            int below = 0;
            double inv_below = 0.0;
            for (int i = r - 1; i > 0; i -= i & (-i)) {
                if (stamp[i] == b) { below += cnt[i]; inv_below += inv[i]; }
            }
            S += (below + 1.0) / r + (total_inv - inv_below);
            for (int i = r; i <= N; i += i & (-i)) {
                if (stamp[i] != b) { stamp[i] = b; cnt[i] = 0; inv[i] = 0.0; }
                cnt[i] += 1;
                inv[i] += 1.0 / r;
            }
            total_inv += 1.0 / r;
            int col = col_of[n];
            if (col > 0) out(b, col - 1) = S / (n + 1);
        }
        for (int n = max_n - 1; n >= 0; n--) std::swap(perm[n], perm[swaps[n]]);
    }
    return out;
}

// Adds to counts[i] the number of null values (column size_col[i] of
// null_ap) at or above ap_obs[i], and keeps in each column of top the
// n_keep largest values seen so far (ascending). counts and top are updated
// in place.
// Counting: the observed values of each size are sorted; each null value v
// is located among them (number of observed values <= v), and a suffix sum
// gives, for every observed value, the number of null values at or above
// it. Tail: only null values above the smallest kept value enter the merge
// with the kept values.
// [[Rcpp::export(name = ".pr_count_chunk", rng = false)]]
void pr_count_chunk(NumericMatrix null_ap, NumericVector ap_obs, IntegerVector size_col,
                    NumericVector counts, NumericMatrix top) {
    int B = null_ap.nrow(), n_sizes = null_ap.ncol(), n_keep = top.nrow();
    std::vector< std::vector<int> > by_size(n_sizes);
    for (int i = 0; i < ap_obs.size(); i++) by_size[size_col[i] - 1].push_back(i);
    std::vector<double> obs_sorted, cand, kept, hits;
    std::vector<int> order;
    for (int k = 0; k < n_sizes; k++) {
        const double* col = &null_ap(0, k);
        const std::vector<int>& idx = by_size[k];
        int m = idx.size();
        if (m > 0) {
            order.assign(idx.begin(), idx.end());
            std::sort(order.begin(), order.end(),
                      [&](int a, int b) { return ap_obs[a] < ap_obs[b]; });
            obs_sorted.resize(m);
            for (int j = 0; j < m; j++) obs_sorted[j] = ap_obs[order[j]];
            hits.assign(m + 1, 0.0);
            for (int b = 0; b < B; b++) {
                int pos = std::upper_bound(obs_sorted.begin(), obs_sorted.end(), col[b]) - obs_sorted.begin();
                hits[pos] += 1.0;
            }
            double running = 0.0;
            for (int j = m; j >= 1; j--) {
                running += hits[j];
                counts[order[j - 1]] += running;
            }
        }
        double threshold = top(0, k);
        cand.clear();
        for (int b = 0; b < B; b++) if (col[b] > threshold) cand.push_back(col[b]);
        if (!cand.empty()) {
            kept.assign(&top(0, k), &top(0, k) + n_keep);
            kept.insert(kept.end(), cand.begin(), cand.end());
            std::nth_element(kept.begin(), kept.end() - n_keep, kept.end());
            std::sort(kept.end() - n_keep, kept.end());
            std::copy(kept.end() - n_keep, kept.end(), &top(0, k));
        }
    }
}
