#include "../kontsevich_graph_series.hpp"
#include <ginac/ginac.h>
#include <iostream>
using namespace std;
using namespace GiNaC;

size_t order = 2;

int main()
{
    // Compute relevant primes
    map< size_t, set<KontsevichGraph> > primes;
    for (size_t n = 0; n <= order; ++n)
    {
        primes[n] = KontsevichGraph::graphs(n, 2, true, true,
                        [](KontsevichGraph g) -> bool
                        {
                            return g.positive_differential_order() && g.is_prime();
                        });
    }
    // Make a table of (symbolic) weights
    map<KontsevichGraph, symbol> weights;
    size_t weight_count = 0;
    for (size_t n = 0; n <= order; ++n)
    {
        for (KontsevichGraph g : primes[n])
        {
            weights[g] = symbol("w_" + to_string(weight_count++));
        }
    }
    KontsevichGraphSeries<symbol> star_product;
}
