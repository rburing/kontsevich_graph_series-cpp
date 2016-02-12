#include "../kontsevich_graph_series.hpp"
#include "../util/partitions.hpp"
#include "../util/factorial.hpp"
#include <ginac/ginac.h>
#include <iostream>
#include <vector>
using namespace std;
using namespace GiNaC;

size_t order = 3;

int main()
{
    cout << "Listing relevant graphs (prime, positive differential order, modulo mirror images) at order";
    map< size_t, set<KontsevichGraph> > relevants;
    for (size_t n = 1; n <= order; ++n)
    {
        cout << " " << n << ",";
        cout.flush();
        relevants[n] = KontsevichGraph::graphs(n, 2, true, true,
                        [](KontsevichGraph g) -> bool
                        {
                            return g.positive_differential_order() && g.is_prime();
                        });
    }
    map< size_t, vector<KontsevichGraph> > primes;
    cout << " done." << endl;
    cout << "Making a table of primes and symbolic weights...";
    map<KontsevichGraph, ex> weights;
    size_t weight_count = 0;
    for (size_t n = 1; n <= order; ++n)
    {
        primes[n].reserve(2*relevants[n].size());
        for (KontsevichGraph g : relevants[n])
        {
            weights[g] = symbol("w_" + to_string(weight_count++));
            primes[n].push_back(g);
            KontsevichGraph mirror = g.mirror_image();
            if (g.abs() != mirror.abs())
            {
                weights[mirror] = weights[g];
                if (n % 2 == 1)
                    weights[mirror] *= -1;
                primes[n].push_back(mirror);
            }
        }
    }
    cout << " done." << endl;
    cout << "Constructing star product at order";
    KontsevichGraphSeries<ex> star_product;
    star_product[0] = KontsevichGraphSum<ex> { { 1, KontsevichGraph(0, 2, {}) } };
    for (size_t n = 1; n <= order; ++n)
    {
        cout << " " << n << ",";
        cout.flush();
        ex major_coeff = 1/(ex)factorial(n);
        Partitions partitions(n);
        for (auto partition = partitions.begin(); partition != partitions.end(); ++partition)
        {
            std::vector<size_t> prime_sizes((*partition).size());
            for (size_t i = 0; i != (*partition).size(); ++i)
                prime_sizes[i] = primes[(*partition)[i]].size();
            CartesianProduct decompositions(prime_sizes);
            for (auto decomposition = decompositions.begin(); decomposition != decompositions.end(); ++decomposition)
            {
                KontsevichGraph composite(0, 2, {});
                ex coeff = major_coeff;
                for (size_t i = 0; i != (*partition).size(); ++i)
                {
                    composite *= primes[(*partition)[i]][(*decomposition)[i]];
                    coeff *= weights[primes[(*partition)[i]][(*decomposition)[i]]];
                }
                coeff *= composite.multiplicity();
                star_product[n] += KontsevichGraphSum<ex>({ { coeff, composite} });
            }
        }
        // TODO: n is not yet considered a partition of n, so treat it separately for now:
        for (KontsevichGraph prime : primes[n])
        {
            star_product[n] += KontsevichGraphSum<ex>({ { major_coeff * weights[prime] * prime.multiplicity(), prime } });
        }
    }
    cout << " done." << endl;
    cout << endl;
    KontsevichGraphSeries<ex> arg = { { 0, { { 1, KontsevichGraph(0, 1, {}) } }} };
    cout << star_product << "\n";
    cout << "Computing associator...";
    cout.flush();
    KontsevichGraphSeries<ex> assoc = star_product({ arg, star_product }) - star_product({ star_product, arg });
    cout << endl;
    cout << "Reducing associator...";
    cout.flush();
    assoc.reduce();
    cout << endl;
    cout << "Number of terms in associator:\n";
    for (size_t n = 0; n <= order; ++n)
    {
        cout << "h^" << n << ":\n";
        cout << assoc[n].size() << " total\n";
        for (std::vector<size_t> indegrees : assoc[n].in_degrees())
        {
            for (size_t j = 0; j != indegrees.size(); ++j)
                cout << indegrees[j] << " ";
            cout << "\n";
        }
    }
}
