#include <Rcpp.h>
#include <algorithm>
using namespace Rcpp;

// Information content of the most informative common ancestor for every
// pair of terms. The ancestors of term k are the entries
// i[p[k]] .. i[p[k + 1] - 1] (ascending term indices, 0-based) of a
// compressed sparse column matrix with one column per term; ic holds the
// information content of every term of the ontology.
// [[Rcpp::export(name = ".mica_ic", rng = false)]]
NumericMatrix mica_ic(IntegerVector p, IntegerVector i, NumericVector ic) {
    int n = p.size() - 1;
    NumericMatrix out(n, n);
    for (int a = 0; a < n; a++) {
        for (int b = a; b < n; b++) {
            int x = p[a], y = p[b];
            const int xe = p[a + 1], ye = p[b + 1];
            double best = 0.0;
            while (x < xe && y < ye) {
                if (i[x] < i[y]) x++;
                else if (i[x] > i[y]) y++;
                else {
                    best = std::max(best, ic[i[x]]);
                    x++;
                    y++;
                }
            }
            out(a, b) = best;
            out(b, a) = best;
        }
    }
    return out;
}
